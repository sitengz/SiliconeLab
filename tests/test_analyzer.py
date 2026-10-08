#!/usr/bin/env python3
"""Scientific contract tests for static analysis and external Z1 integration."""
import argparse
import csv
import importlib.util
import json
import math
from pathlib import Path
import subprocess
import sys
import tempfile
import unittest

ARGS = None


class AnalyzerTests(unittest.TestCase):
    def setUp(self):
        self.tmp = tempfile.TemporaryDirectory()
        self.root = Path(self.tmp.name)

    def tearDown(self):
        self.tmp.cleanup()

    def generate(self, system, extra=""):
        cfg = self.root / "model.conf"
        common = {"oil": "chain_count = 4 3\nmps_percent = 0\n", "elastomer": "strand_length = 8\nstrand_count = 4\ncrosslinker_length = 8\nmoderator_count = 0\n",
                  "coating": "model = v35\nn1 = 8\nm1 = 4\nn2 = 8\nm4 = 6\nfunctionality = 4\noil = copolymer\noil_length = 5\nmps_percent = 50\noil_wt = 30\nmps_distribution = balanced\n"}[system]
        cfg.write_text(f"system = {system}\noutput_dir = run\n" + common + extra)
        subprocess.run([ARGS.generator, "--config", str(cfg)], check=True, capture_output=True)
        run = self.root / "run"
        return run / "data.run", run / "run.info"

    def analyze(self, data, info, name="analysis", extra=(), success=True):
        out = self.root / name
        p = subprocess.run([sys.executable, ARGS.driver, str(data), str(info),
                            "--backend-dir", ARGS.backend_dir, "--output-dir", str(out), *extra],
                           capture_output=True, text=True)
        if success:
            self.assertEqual(p.returncode, 0, p.stdout + p.stderr + ((out / "analysis.log").read_text() if (out / "analysis.log").exists() else ""))
        return p, out

    def read(self, path):
        with path.open() as f:
            return list(csv.DictReader(f, delimiter="\t"))

    def test_pure_oil_network_not_applicable_and_export(self):
        data, info = self.generate("oil")
        _, out = self.analyze(data, info, extra=("--z1-mode", "export"))
        report = json.loads((out / "analysis.json").read_text())
        for i in (3, 4, 5, 6, 7, 10):
            self.assertEqual(report["analyses"][i-1]["status"], "not_applicable")
        mapping = json.loads((out / "z1/chain_mapping.json").read_text())
        self.assertEqual(len(mapping["chains"]), 3)
        self.assertEqual(len(self.read(out / "oil_chain_properties.tsv")), 3)

    def test_known_oil_geometry_across_periodic_boundary(self):
        data, info = self.generate("oil")
        text = data.read_text()
        # Preserve complete topology, but replace first chain with known wrapped positions.
        import re
        match = re.search(r"([\d.eE+-]+) ([\d.eE+-]+) xlo xhi", text)
        lo, hi = map(float, match.groups())
        length = hi-lo
        lines=text.splitlines()
        section=""
        for n,line in enumerate(lines):
            if line.startswith("Atoms"): section="atoms";continue
            if line.startswith(("Bonds", "Angles", "Dihedrals", "Velocities")): section="other"
            fields=line.split()
            if section=="atoms" and len(fields)>=7 and fields[0].isdigit() and int(fields[1])==1:
                rank=int(fields[0])-1
                x=hi-1+rank
                fields[4]=str(lo+(x-lo)%length);fields[5]="0";fields[6]="0"
                lines[n]=" ".join(fields)
        data.write_text("\n".join(lines)+"\n")
        _, out = self.analyze(data,info,extra=("--z1-mode","export"))
        row=self.read(out / "oil_chain_properties.tsv")[0]
        self.assertAlmostEqual(float(row["Lc_A"]),3)
        self.assertAlmostEqual(float(row["Ree_A"]),3)
        self.assertAlmostEqual(float(row["Rg_A"]),math.sqrt(1.25))
        self.assertAlmostEqual(float(row["anisotropy"]),1)

    def test_balanced_mps_counts_and_backbone_selection(self):
        data,info=self.generate("coating")
        _,out=self.analyze(data,info,extra=("--z1-mode","export"))
        rows=self.read(out / "oil_chain_properties.tsv")
        metadata=json.loads(info.read_text())
        self.assertTrue(rows)
        self.assertEqual(sum(int(r["beads"]) for r in rows),metadata["components"]["silicone_oil_filler"]["beads"])
        self.assertTrue(all(int(r["backbone_beads"])==5 for r in rows))
        self.assertTrue(all(int(r["beads"])>int(r["backbone_beads"]) for r in rows))

    def test_mismatched_metadata_and_preserved_input(self):
        data,info=self.generate("oil")
        before=data.read_bytes()
        d=json.loads(info.read_text());d["topology_counts"]["atoms"]+=1
        info.write_text(json.dumps(d))
        p,out=self.analyze(data,info,extra=("--z1-mode","export"),success=False)
        self.assertNotEqual(p.returncode,0)
        self.assertEqual(data.read_bytes(),before)

    def test_triclinic_box_rejected(self):
        data,info=self.generate("oil")
        data.write_text(data.read_text().replace("Masses", "1 0 0 xy xz yz\n\nMasses", 1))
        p,out=self.analyze(data,info,extra=("--z1-mode","export"),success=False)
        self.assertNotEqual(p.returncode,0)
        self.assertIn("triclinic",(out/"analysis.log").read_text())

    def test_elastomer_parity_with_native_backend(self):
        data,info=self.generate("elastomer")
        # A fully reacted two-junction/four-strand graph tests nontrivial reduction,
        # parallel edges and modulus tables without launching molecular dynamics.
        metadata=json.loads(info.read_text())
        strands=metadata["components"]["strands"]
        crosslinkers=metadata["components"]["crosslinkers"]
        ends=[[m*strands["N"]+s for s in metadata["strand"]["reactive_bead_sites"]]
              for m in range(strands["M"])]
        offset=strands["beads"]
        sites=[[offset+m*crosslinkers["N"]+s
                for s in metadata["crosslinker"]["reactive_bead_sites"]]
               for m in range(crosslinkers["M"])]
        pairs=[(ends[m][j],sites[j][m]) for m in range(strands["M"]) for j in (0,1)]
        reacted={a for pair in pairs for a in pair}
        lines=data.read_text().splitlines();section="";bond_count=metadata["topology"]["bonds"]
        for i,line in enumerate(lines):
            fields=line.split()
            if len(fields)==2 and fields[1]=="bonds":lines[i]=f"{bond_count+len(pairs)} bonds"
            if line.startswith("Atoms"):section="atoms";continue
            if line.startswith(("Bonds","Angles","Dihedrals","Velocities")):section="other"
            if section=="atoms" and len(fields)>=7 and fields[0].isdigit() and int(fields[0]) in reacted:
                fields[2]="1";lines[i]=" ".join(fields)
        start=next(i for i,line in enumerate(lines) if line.startswith("Bonds"))+1
        stop=next((i for i in range(start,len(lines)) if lines[i].startswith(("Angles","Dihedrals","Velocities"))),len(lines))
        lines[stop:stop]=[f"{bond_count+i+1} 2 {a} {b}" for i,(a,b) in enumerate(pairs)]+[""]
        data.write_text("\n".join(lines)+"\n")
        _,out=self.analyze(data,info,extra=("--z1-mode","export"))
        report=json.loads((out/"analysis.json").read_text())
        self.assertEqual(report["chemical_network"]["overall_conversion"],1)
        self.assertEqual(report["analyses"][10]["chains"],4)
        native=self.root / "native"
        subprocess.run([ARGS.reference_basic,str(data),str(info),"--output-dir",str(native)],check=True,capture_output=True)
        for name in ("strand_properties.run.tsv","network_statistics.run.tsv","strand_statistics.run.tsv","junction_properties.run.tsv"):
            self.assertEqual((out/name).read_bytes(),(native/name).read_bytes(),name)

    def test_z1_zero_exit_without_outputs_is_failure(self):
        data,info=self.generate("oil")
        fake=self.root / "fake_z1"
        fake.write_text("#!/bin/sh\nexit 0\n");fake.chmod(0o755)
        p,out=self.analyze(data,info,extra=("--z1-command",str(fake)),success=False)
        self.assertNotEqual(p.returncode,0)
        self.assertEqual(json.loads((out/"analysis.json").read_text())["status"],"failed")

    def test_validated_z1_fixture_and_mapping(self):
        data,info=self.generate("oil")
        fake=self.root / "fake_z1"
        fake.write_text('''#!/usr/bin/env python3
from pathlib import Path
lines=Path('config.Z1').read_text().splitlines()
n=int(lines[0]);lengths=list(map(int,lines[2].split()));box=lines[1]
import math
cursor=3;lp=[];ree=[]
with Path('Z1+SP.dat').open('w') as f:
 f.write(str(n)+'\\n'+box+'\\n')
 for size in lengths:
  points=[list(map(float,l.split())) for l in lines[cursor:cursor+size]];cursor+=size
  f.write('2\\n')
  for p,s in ((points[0],1),(points[-1],size)): f.write(' '.join(map(str,p))+f' {s} 0\\n')
  lp.append(math.dist(points[0],points[-1]));ree.append(lp[-1]**2)
for name,values in [('Lpp_values.dat',lp),('Z_values.dat',[0]*n),('N_values.dat',lengths)]: Path(name).write_text(' '.join(map(str,values))+'\\n')
Path('Z1+summary.dat').write_text(' '.join(map(str,[0,n,sum(lengths)/n,sum(ree)/n,sum(lp)/n,0,0,0,0,float('inf'),float('nan'),float('inf'),float('nan'),1,1]))+'\\n')
''');fake.chmod(0o755)
        _,out=self.analyze(data,info,extra=("--z1-command",str(fake)))
        report=json.loads((out/"analysis.json").read_text())
        self.assertEqual(report["analyses"][11]["status"],"complete")
        self.assertEqual(report["primitive_paths"]["groups"]["oil"]["chains"],3)

    def test_existing_output_is_not_overwritten(self):
        data,info=self.generate("oil")
        self.analyze(data,info,extra=("--z1-mode","export"))
        p,out=self.analyze(data,info,extra=("--z1-mode","export"),success=False)
        self.assertNotEqual(p.returncode,0)


if __name__=="__main__":
    p=argparse.ArgumentParser()
    for name in ("generator","driver","backend-dir","reference-basic"): p.add_argument("--"+name,required=True)
    ARGS,rest=p.parse_known_args()
    unittest.main(argv=[sys.argv[0],*rest])
