"""Composite contracts across native matrix backends and independent NP pairs."""
import argparse
from collections import defaultdict
import csv
import json
import math
from pathlib import Path
import subprocess
import tempfile
import unittest

p=argparse.ArgumentParser()
p.add_argument('--generator', type=Path, required=True)
ARGS,REST=p.parse_known_args()
GEN=ARGS.generator.resolve()

def structure(path):
    text=path.read_text()
    box=[]
    for axis in 'xyz':
        r=next(x.split() for x in text.splitlines() if x.endswith(f'{axis}lo {axis}hi'))
        box.append(float(r[1])-float(r[0]))
    masses={int(r[0]):float(r[1]) for line in text.split('Masses\n\n')[1].split('\n\nAtoms')[0].splitlines() if (r:=line.split())}
    atoms={}
    for line in text.split('Atoms # full\n\n')[1].split('\nBonds\n')[0].splitlines():
        if not (r:=line.split()): continue
        ident,mol,typ=map(int,r[:3]); q=float(r[3])
        xyz=[float(r[4+j])+(int(r[7+j])*box[j] if len(r)==10 else 0) for j in range(3)]
        atoms[ident]=(mol,typ,q,xyz)
    def section(name,next_name):
        if '\n'+name+'\n' not in text: return []
        block=text.split('\n'+name+'\n')[1].split('\n'+next_name+'\n')[0]
        return [list(map(int,r.split())) for r in block.splitlines() if r.strip()]
    return text,box,masses,atoms,section('Bonds','Angles'),section('Angles','Dihedrals')

class CompositeContracts(unittest.TestCase):
    def setUp(self):
        self.temp=tempfile.TemporaryDirectory()
        self.addCleanup(self.temp.cleanup)
        self.root=Path(self.temp.name)

    def generate(self,base='v22',extra='',name='case',system='composite',ok=True,radius=10,count=2):
        if base in ('v22','v35'):
            native='n1=24\nm1=24\nn2=16\nfunctionality=4\nm4=0\n'
            select=f'model={base}\n' if system=='coating' else f'base={base}\n'
        else:
            native='strand_length=16\nstrand_count=24\ncrosslinker_length=16\nfunctionality=4\nmoderator_count=0\n'
            if 'strand_topology=grafted' in extra: native=native.replace('strand_length=16\n','')
            select='' if system=='elastomer' else f'base={base}\n'
        conf=self.root/(name+'.conf')
        conf.write_text(f'system={system}\n'+select+f'output_dir={name}\n'+native+
                        f'density=0.1\ntarget_density=0.8\nseed=5489\nnp_count={count}\nnp_radius={radius}\n'+extra)
        result=subprocess.run([str(GEN),'--config',str(conf)],capture_output=True,text=True,timeout=60)
        self.assertEqual(result.returncode==0,ok,result.stdout+result.stderr)
        return self.root/name,result

    def audit(self,run):
        meta=json.loads((run/'composite.json').read_text())
        native,bbox,bmass,b,bbonds,bangles=structure(run/('precursor/data.'+run.name))
        actual,box,masses,a,bonds,angles=structure(run/('data.'+run.name))
        self.assertEqual(native.split('\nBonds\n')[1],actual.split('\nBonds\n')[1])
        self.assertEqual({i:a[i][:3] for i in b},{i:v[:3] for i,v in b.items()})
        self.assertEqual({i:masses[i] for i in bmass},bmass)
        np_type=max(bmass)+1
        self.assertEqual(meta['np_atom_type'],np_type)
        groups=defaultdict(list)
        for i,v in a.items():
            if i in b: continue
            self.assertEqual(v[1:3],(np_type,0.))
            groups[v[0]].append(v[3])
        self.assertEqual(len(groups),meta['np_count'])
        self.assertEqual(set(map(len,groups.values())),{meta['beads_per_np']})
        matrix_mass=sum(bmass[v[1]] for v in b.values())
        np_mass=masses[np_type]*sum(map(len,groups.values()))
        self.assertAlmostEqual(meta['np_weight_percent'],100*np_mass/(matrix_mass+np_mass))
        self.assertAlmostEqual(box[0]**3/bbox[0]**3,(matrix_mass+np_mass)/matrix_mass)
        with (run/'np-template.csv').open() as f: template=list(csv.DictReader(f))
        self.assertEqual(len(template),meta['beads_per_np'])
        pts=[[float(row[k]) for k in ('x_A','y_A','z_A')] for row in template]
        self.assertAlmostEqual(min(math.dist(x,y) for i,x in enumerate(pts) for y in pts[i+1:]),meta['fcc_nearest_neighbor_A'],places=10)
        self.assertTrue(all(sum(x*x for x in q)<=meta['np_radius_bead_centers_A']**2+1e-9 for q in pts))
        def vec(points,i,j,periodic=None):
            delta=[points[i][3][k]-points[j][3][k] for k in range(3)]
            return [d-periodic[k]*round(d/periodic[k]) for k,d in enumerate(delta)] if periodic else delta
        def norm(v): return math.sqrt(sum(x*x for x in v))
        for _,_,i,j in bbonds:
            self.assertAlmostEqual(norm(vec(b,i,j,bbox)),norm(vec(a,i,j)),places=7)
        # Bond and angle invariants catch incorrectly unwrapped branches/rings.
        for _,_,i,j,k in bangles:
            def cosine(points,periodic):
                u,v=vec(points,i,j,periodic),vec(points,k,j,periodic)
                return sum(x*y for x,y in zip(u,v))/(norm(u)*norm(v))
            self.assertAlmostEqual(cosine(b,bbox),cosine(a,None),places=7)
        for i in b:
            for c in meta['np_centers_A']:
                d=math.sqrt(sum((a[i][3][j]-c[j]-box[j]*round((a[i][3][j]-c[j])/box[j]))**2 for j in range(3)))
                self.assertGreaterEqual(d+1e-8,meta['matrix_np_center_clearance_A'])
        before=(run/('precursor/in.'+run.name)).read_text().splitlines()
        after=(run/('in.'+run.name)).read_text().splitlines()
        for cmd in ('bond_coeff','angle_coeff','dihedral_coeff','timestep','minimize'):
            self.assertEqual([x.split() for x in before if x.split()[:1]==[cmd]],
                             [x.split() for x in after if x.split()[:1]==[cmd]])
        precursor=[]
        np_pairs=[]
        hot=True
        for line in after:
            r=line.split()
            if r[:2]==['pair_style','hybrid']: hot=r[2]=='lj/cut'
            if r[:1]!=['pair_coeff']: continue
            values=r[5 if hot else 4:]
            if int(r[2])<np_type: precursor.append(r[:3]+values)
            else: np_pairs.append(list(map(float,values)))
        self.assertEqual(precursor,[x.split() for x in before if x.split()[:1]==['pair_coeff']])
        self.assertEqual(len(np_pairs),2*np_type)
        for phase in (0,1):
            for coeff,pair in zip(np_pairs[phase*np_type:(phase+1)*np_type],meta['np_pairs']):
                self.assertEqual(coeff,[pair['epsilon_kcal_mol'],pair['sigma_A'],pair['cutoff_A']])
                eps,sigma,cut=coeff
                u=lambda r:4*eps*((sigma/r)**12-(sigma/r)**6)-4*eps*((sigma/cut)**12-(sigma/cut)**6)
                self.assertAlmostEqual(u(cut),0)
                if pair['interaction']=='wca':
                    self.assertAlmostEqual(cut,2**(1/6)*sigma)
                    self.assertGreater(u(.95*sigma),0)
                else:
                    self.assertGreater(cut,2**(1/6)*sigma)
                    self.assertLess(u(min(1.5*sigma,(cut+2**(1/6)*sigma)/2)),0)
        record=json.loads((run/'siliconelab.json').read_text())
        self.assertEqual(record['base'],meta['base'])
        self.assertEqual(record['series'],'composite')
        self.assertEqual(record['matrix_system'],meta['matrix_system'])
        self.assertEqual(set(record['generated_files']),{p.name for p in run.iterdir() if p.is_file()})
        self.assertFalse(any(x.split()[:4] in (['fix','integrate','all','nvt'],['fix','integrate','all','npt']) for x in after))
        self.assertEqual((run/'nanoparticles.json').read_bytes(),(run/'composite.json').read_bytes())
        return meta,after

    def test_all_bases_and_oil_filled_coatings(self):
        for base in ('v22','v35','elastomer'):
            with self.subTest(base=base):
                run,_=self.generate(base,name=base)
                meta,_=self.audit(run)
                self.assertEqual(meta['beads_per_np'],13)
        for oil in ('pdms','pmps','copolymer'):
            with self.subTest(oil=oil):
                extra=f'oil={oil}\noil_length=4\noil_wt=10\n'
                if oil=='copolymer':extra+='mps_percent=50\n'
                run,_=self.generate('v35',extra=extra,name=oil)
                self.audit(run)

    def test_elastomer_architectures_and_conversion(self):
        for topology,extra in (
            ('ring','strand_functionality=4\n'),
            ('star','strand_length=8\nstrand_arm_count=4\n'),
            ('grafted','backbone_length=12\nside_chain_length=4\ngraft_spacing=2\ngraft_functional_fraction=50\n')):
            with self.subTest(topology=topology):
                # Grafted backends forbid even an otherwise unused strand_length.
                run,_=self.generate('elastomer',extra=f'strand_topology={topology}\n'+extra,name=topology)
                meta,after=self.audit(run)
                self.assertIn('fix np_minimize_hold nanoparticles setforce 0 0 0',after)
                self.assertTrue(any(x.startswith('velocity precursor create ') for x in after))
        run,_=self.generate('elastomer',extra='target_conversion=80\nfiller_length=4\nfiller_wt=10\n',name='controlled')
        _,after=self.audit(run)
        self.assertEqual(sum(x.startswith('run 1 start ') for x in after),4)
        self.assertTrue(all('$(step+1:%.0f)' in x for x in after if x.startswith('run 1 start ')))
        self.assertTrue(any(x.split()[:4]==['fix','conversion_halt','all','halt'] for x in after))

    def test_independent_interactions_size_and_masses(self):
        run,_=self.generate('elastomer',radius=15,count=1,name='size')
        meta,_=self.audit(run)
        self.assertEqual(meta['beads_per_np'],55)
        self.assertEqual(meta['np_count'],1)
        for i,modes in enumerate((('wca','lj'),('lj','wca'))):
            with self.subTest(modes=modes):
                extra=(f'np_np_interaction={modes[0]}\nnp_matrix_interaction={modes[1]}\n'
                       'np_sigma=7\nnp_epsilon=3\nnp_np_epsilon=2\n'
                       'np_matrix_epsilon=0.8\nnp_matrix_sigma=6.2\nnp_bead_mass=148\nmass=100\n')
                if modes[0]=='lj':extra+='np_np_cutoff=20\n'
                if modes[1]=='lj':extra+='np_matrix_cutoff=16\n'
                run,_=self.generate('elastomer',radius=15,count=3,extra=extra,name=f'pairs-{i}')
                meta,_=self.audit(run)
                self.assertEqual(meta['np_np_interaction'],modes[0])
                self.assertEqual(meta['np_matrix_interaction'],modes[1])
                self.assertEqual(meta['np_pairs'][-1]['epsilon_kcal_mol'],2)
                self.assertEqual(meta['np_pairs'][0]['epsilon_kcal_mol'],.8)
                self.assertEqual(meta['np_pairs'][0]['sigma_A'],6.2)
                self.assertEqual(meta['mass_g_mol_per_bead'],148)

    def test_determinism_and_matrix_separation(self):
        first,_=self.generate(name='first',extra='np_spacing=8\nnp_seed=31415\n')
        second,_=self.generate(name='second',extra='np_spacing=8\nnp_seed=31415\n')
        _,box,_,atoms,_,_=structure(first/'data.first')
        self.assertEqual(atoms,structure(second/'data.second')[3])
        matrix=list(atoms.values())[:-26]
        # Independent exhaustive check on a small fixture; no shared grid logic.
        closest=min(math.sqrt(sum((a[3][k]-b[3][k]-box[k]*round((a[3][k]-b[3][k])/box[k]))**2 for k in range(3)))
                    for i,a in enumerate(matrix) for b in matrix[i+1:] if a[0]!=b[0])
        self.assertGreaterEqual(closest,4.5-1e-9)
        self.audit(first)

    def test_interaction_only_changes_preserve_initial_structure(self):
        wca,_=self.generate(name='repulsive')
        lj,_=self.generate(name='attractive',extra='np_interaction=lj\n')
        self.assertEqual(structure(wca/'data.repulsive')[3],structure(lj/'data.attractive')[3])
        self.audit(lj)

    def test_capabilities_and_legacy_aliases(self):
        result=subprocess.run([str(GEN),'--capabilities'],capture_output=True,text=True,check=True)
        self.assertEqual(json.loads(result.stdout)['composite']['bases'],['v22','v35','elastomer'])
        capability=json.loads(result.stdout)['composite']
        self.assertTrue(capability['analysis_supported'])
        self.assertFalse(capability['np_obstacle_primitive_paths_supported'])
        self.assertFalse(capability['dynamics_analysis_supported'])
        result=subprocess.run([str(GEN),'--help','composite'],capture_output=True,text=True,check=True)
        self.assertIn('np_matrix_interaction',result.stdout)
        for system,base in (('coating','v22'),('coating','v35'),('elastomer','elastomer')):
            run,_=self.generate(base,system=system,name='alias-'+base)
            self.audit(run)

    def test_invalid_options_and_no_overwrite(self):
        for base in ('v99','oil'):
            run,_=self.generate(base=base,name='bad-base-'+base,ok=False)
            self.assertFalse(run.exists())
        for i,extra in enumerate(('np_count=3','np_interaction=unknown','np_matrix_interaction=unknown',
                                  'np_np_cutoff=20','np_interaction=lj\nnp_np_cutoff=1',
                                  'np_bead_mass=nan','np_sigma=-1','np_representation=giant',
                                  'np_matrix_sigma=1000','np_spacing=100','thickness=30')):
            run,_=self.generate(extra=extra,name=f'bad-{i}',ok=False)
            self.assertFalse(run.exists())
        run,_=self.generate()
        original=(run/'data.case').read_bytes()
        self.generate(ok=False)
        self.assertEqual((run/'data.case').read_bytes(),original)

if __name__=='__main__':
    unittest.main(argv=[__file__,*REST])
