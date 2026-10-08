#!/usr/bin/env python3
"""Scientific stage durations, coefficient selection and dependency/output guards."""
import importlib.util
import json
from pathlib import Path
import tempfile
import unittest
ROOT=Path(__file__).resolve().parents[1]
spec=importlib.util.spec_from_file_location('workflow',ROOT/'scripts/film_workflow.py')
w=importlib.util.module_from_spec(spec);spec.loader.exec_module(w)

class FilmContracts(unittest.TestCase):
    def test_force_field_selection(self):
        with tempfile.TemporaryDirectory() as d:
            p=Path(d)/'in.native';p.write_text('units real\natom_style full\nbond_style harmonic\nangle_style harmonic\npair_style lj/cut 7\nread_data data.init extra/bond/per/atom 4\nbond_coeff 1 115 2.8\nbond_coeff 2 115 7.2\npair_coeff 1 1 0.5 6.6 7.4\npair_style lj/gromacs 12 15\nbond_coeff 2 115 2.8\npair_coeff 1 1 1.01 6.44\n')
            text=w.cold_header(p,'../data.final')
            self.assertIn('pair_style lj/gromacs 12 15',text);self.assertNotIn('lj/cut',text)
            self.assertIn('bond_coeff 2 115 2.8',text);self.assertNotIn('115 7.2',text)
            self.assertLess(text.index('read_data'),text.index('bond_coeff'))
    def test_sampling_and_independent_origins(self):
        with tempfile.TemporaryDirectory() as d:
            p=Path(d)/'in';w.build_dynamics(p,'# native header\n');t=p.read_text()
            self.assertIn('custom 20 dump.debye',t);self.assertIn('run 20000',t)
            self.assertIn('run 4000000',t);self.assertEqual(t.count('reset_timestep 0 time 0.0'),2)
            w.build_tensile(p,'# native header\n','x',100);t=p.read_text()
            self.assertIn('y 1 1 500',t);self.assertIn('x erate 2e-9',t);self.assertIn('run 20000000',t)
    def test_partial_parent_blocks_children(self):
        with tempfile.TemporaryDirectory() as d:
            old=w.ROOT;w.ROOT=Path(d)
            try:
                with self.assertRaises(RuntimeError):w.check_stage({'directory':'.','required_outputs':[],'final_data':'missing','expected_atoms':1})
            finally:w.ROOT=old
    def test_matched_compositions(self):
        plan=json.loads((ROOT/'workflows/films.json').read_text());self.assertEqual(len(plan['cases']),4)
        for c in plan['cases']:
            film=(ROOT/c['config']).read_text()
            bulk=(ROOT/'examples/full-size'/Path(c['config']).name).read_text()
            scientific=lambda t:[x.strip() for x in t.splitlines() if x.strip() and not x.startswith('#') and not x.startswith(('output_dir','thickness','film_padding'))]
            self.assertEqual(scientific(film),scientific(bulk))
if __name__=='__main__':unittest.main()
