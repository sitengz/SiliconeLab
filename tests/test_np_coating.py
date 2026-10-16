"""NP extension contracts: independent topology, geometry and force-field checks."""
import argparse
from collections import Counter, defaultdict
import importlib.util
import json
import math
from pathlib import Path
import subprocess
import tempfile
import unittest

p=argparse.ArgumentParser()
p.add_argument('--generator',type=Path,required=True)
ARGS,REST=p.parse_known_args()
GEN=ARGS.generator.resolve()

def atoms(path):
    text=path.read_text()
    block=text.split('Atoms # full\n\n')[1].split('\nBonds\n')[0]
    return [line.split() for line in block.splitlines() if line.strip()]

class NPContracts(unittest.TestCase):
    def setUp(self):
        self.temp=tempfile.TemporaryDirectory()
        self.addCleanup(self.temp.cleanup)
        self.root=Path(self.temp.name)

    def generate(self,extra='',name='np',ok=True):
        conf=self.root/(name+'.conf')
        conf.write_text('system=coating\nmodel=v22\nn1=16\nm1=24\nn2=16\n'
                        'functionality=8\nm4=0\ndensity=0.1\ntarget_density=0.8\n'
                        f'output_dir={name}\n'+extra)
        result=subprocess.run([str(GEN),'--config',str(conf)],capture_output=True,text=True,timeout=60)
        self.assertEqual(result.returncode==0,ok,result.stdout+result.stderr)
        return self.root/name

    def prepared(self,name='np'):
        d=2**(1/6)*6.44584366
        return self.generate(f'np_count=2\nnp_radius={d:.17g}\nnp_seed=20261014\n',name)

    def test_geometry_topology_and_mass(self):
        run=self.prepared()
        base=run/'precursor/data.np'
        final=run/'data.np'
        b=atoms(base); a=atoms(final)
        info=json.loads((run/'np.info').read_text())
        self.assertEqual(len(b),480)
        self.assertEqual(len(a),506)
        self.assertEqual(info['beads_per_np'],13)
        self.assertEqual(info['np_count'],2)
        self.assertAlmostEqual(info['np_weight_percent'],100*26/506)
        self.assertEqual(final.read_text().split('\nBonds\n')[1],base.read_text().split('\nBonds\n')[1])
        self.assertEqual([[float(v) for v in r[:4]] for r in a[:480]],
                         [[float(v) for v in r[:4]] for r in b])
        np_groups=defaultdict(list)
        L=info['box_length_A']
        for r in a[480:]:
            self.assertEqual(int(r[2]),6)
            self.assertEqual(float(r[3]),0)
            np_groups[int(r[1])].append(tuple(float(x) for x in r[4:7]))
        self.assertEqual(sorted(map(len,np_groups.values())),[13,13])
        d=2**(1/6)*6.44584366
        for points in np_groups.values():
            self.assertAlmostEqual(min(math.dist(x,y) for i,x in enumerate(points) for y in points[i+1:]),d,places=10)
        def xyz(row,box):
            return [float(row[4+j])+(int(row[7+j])*box if len(row)>7 else 0) for j in range(3)]
        bounds=next(line.split() for line in base.read_text().splitlines() if line.endswith('xlo xhi'))
        oldL=float(bounds[1])-float(bounds[0])
        # All native bond lengths survive randomized molecule packing.
        bonds=base.read_text().split('\nBonds\n')[1].split('\nAngles\n')[0]
        for row in bonds.splitlines():
            if not row.strip(): continue
            _,_,i,j=map(int,row.split())
            old=math.dist(xyz(b[i-1],oldL),xyz(b[j-1],oldL))
            new=math.dist(xyz(a[i-1],L),xyz(a[j-1],L))
            self.assertAlmostEqual(old,new,places=7)
        # NP-containing pairs start outside the WCA range, including periodic images.
        for row in a[:480]:
            q=[float(v) for v in row[4:7]]
            for center in info['np_centers_A']:
                dist=math.sqrt(sum((q[j]-center[j]-L*round((q[j]-center[j])/L))**2 for j in range(3)))
                self.assertGreaterEqual(dist,info['np_radius_bead_centers_A']+d-1e-9)
        # LAMMPS Masses needs all six rows without an intervening blank line.
        masses=final.read_text().split('Masses\n\n')[1].split('\n\nAtoms')[0].strip().splitlines()
        self.assertEqual(len(masses),6)
        self.assertEqual(masses[-1].split()[:2],['6','74'])

    def test_native_interactions_and_stages_preserved(self):
        run=self.prepared()
        native=(run/'precursor/in.np').read_text().splitlines()
        result=(run/'in.np').read_text().splitlines()
        for command in ('bond_coeff','angle_coeff','dihedral_coeff','run','timestep'):
            self.assertEqual([x.split() for x in native if x.split()[:1]==[command]],
                             [x.split() for x in result if x.split()[:1]==[command]
                              and not x.startswith('run 1 start ')])
        self.assertEqual(sum(x.startswith('run 1 start ') for x in result),3)
        original=[x.split() for x in native if x.split()[:1]==['pair_coeff']]
        pairs=[x.split() for x in result if x.split()[:1]==['pair_coeff']]
        precursor=[x for x in pairs if int(x[2])<=5]
        restored=[]
        for x in precursor:
            start=5 if x[3]=='lj/cut' else 4
            restored.append(x[:3]+x[start:])
        self.assertEqual(original,restored)
        np_pairs=[x for x in pairs if x[2]=='6']
        self.assertEqual(len(np_pairs),12)
        for x in np_pairs:
            start=5 if x[4]=='2' else 4
            epsilon,sigma,cut=map(float,x[start:])
            self.assertGreater(epsilon,0)
            self.assertAlmostEqual(cut,2**(1/6)*sigma,places=12)
        self.assertIn('neigh_modify exclude molecule/intra nanoparticles',result)
        self.assertIn('compute precursor_temp precursor temp',result)
        self.assertIn('compute composite_press all pressure composite_temp',result)
        self.assertEqual(result.count('fix_modify integrate temp precursor_temp press composite_press'),2)
        self.assertEqual(sum('dilate precursor' in x for x in result),2)
        self.assertTrue(any(x.startswith('fix compress precursor deform ') for x in result))
        self.assertFalse(any(x.split()[0]=='fix' and x.split()[2]=='all' and x.split()[3] in ('nvt','npt')
                             for x in result if x.split()))
        self.assertEqual(sum('rigid/nvt/small' in x for x in result),5)
        self.assertEqual(sum('reinit no' in x for x in result),5)
        active=set()
        for i,line in enumerate(result):
            words=line.split()
            if words[:1]==['unfix']: active.discard(words[1])
            if words[:1]==['fix']: active.add(words[1])
            if line.split()[:4] in (['fix','integrate','precursor','npt'],
                                    ['fix','integrate','precursor','nvt'],
                                    ['fix','msd_integrate','precursor','nvt']):
                self.assertIn('np_'+words[1],active)
                if words[3]=='npt': self.assertIn('rigid/nvt/small',result[i-1])

    def test_determinism_and_no_overwrite(self):
        first=self.prepared('first'); second=self.prepared('second')
        self.assertEqual(atoms(first/'data.first'),atoms(second/'data.second'))
        old=(first/'data.first').read_bytes()
        self.generate('np_count=2\n',name='first',ok=False)
        self.assertEqual(old,(first/'data.first').read_bytes())

    def test_invalid_options(self):
        for i,extra in enumerate(('np_count=0','np_count=513','np_count=2.5','np_radius=30',
                                  'np_count=2\nnp_radius=nan','np_count=2\nnp_radius=1',
                                  'np_count=2\nnp_count=3','np_count=2\nnp_interaction=unknown',
                                  'np_count=2\nthickness=30')):
            folder=self.generate(extra,name=f'bad{i}',ok=False)
            self.assertFalse(folder.exists())

    def test_analysis_recognizes_composite_population(self):
        run=self.prepared()
        spec=importlib.util.spec_from_file_location('driver',Path(__file__).resolve().parents[1]/'src/siliconelab_analyzer.py')
        driver=importlib.util.module_from_spec(spec); spec.loader.exec_module(driver)
        import sys
        sys.path.insert(0,str(Path(__file__).resolve().parents[1]/'src'))
        normalized=driver.normalize(json.loads((run/'np.info').read_text()))
        self.assertEqual(normalized['source_system'],'composite')
        self.assertEqual(normalized['matrix_system'],'coating')
        self.assertLess(normalized['composition']['total_beads'],normalized['composite']['total_beads'])

if __name__=='__main__':
    unittest.main(argv=[__file__,*REST])
