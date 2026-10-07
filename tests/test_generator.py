"""Scientific parity and run-contract checks; no MD jobs are executed."""
import argparse
import hashlib
import json
from pathlib import Path
import subprocess
import tempfile
import unittest

parser = argparse.ArgumentParser()
parser.add_argument("--generator", type=Path, required=True)
parser.add_argument("--reference-dir", type=Path, required=True)
parser.add_argument("--source-dir", type=Path, required=True)
ARGS, UNITTEST_ARGS = parser.parse_known_args()
ARGS.generator = ARGS.generator.resolve()
ARGS.reference_dir = ARGS.reference_dir.resolve()
ARGS.source_dir = ARGS.source_dir.resolve()


def run(executable, *args, cwd):
    return subprocess.run([str(executable), *map(str, args)], cwd=cwd,
                          text=True, capture_output=True, timeout=45)


def example_entries(filename):
    entries = []
    for line in (ARGS.source_dir / "examples" / filename).read_text().splitlines():
        line = line.partition("#")[0].strip()
        if not line:
            continue
        key, value = map(str.strip, line.split("=", 1))
        if key not in ("system", "model", "output_dir"):
            entries.append((key, value))
    return entries


class GeneratorContract(unittest.TestCase):
    def setUp(self):
        self.temp = tempfile.TemporaryDirectory(prefix="SiliconeLab test ")
        self.addCleanup(self.temp.cleanup)
        self.root = Path(self.temp.name)
        self.cwd = self.root / "unrelated working directory"
        self.cwd.mkdir()

    def assert_ok(self, result):
        self.assertEqual(result.returncode, 0, result.stdout + result.stderr)

    def parity(self, system, model, details, name):
        config = self.root / (name + ".conf")
        selector = f"system = {system}\n" + (f"model = {model}\n" if model else "")
        detail_text = "".join(f"{k} = {v}\n" for k, v in details)
        config.write_text(selector + f'output_dir = "unified/{name}/"\n' + detail_text)
        self.assert_ok(run(ARGS.generator, "--config", config, cwd=self.cwd))
        package = self.root / "unified" / name
        ref_parent = self.root / "reference"
        ref_parent.mkdir(exist_ok=True)
        ref_package = ref_parent / name
        native_output = (ref_package if system == "oil" else ref_parent) / ("data." + name)
        native_config = self.root / (name + ".native.conf")
        native_config.write_text(detail_text + f'output = "{native_output.as_posix()}"\n')
        backend_name = system if model is None else "coating_" + model
        suffix = ARGS.generator.suffix
        reference = ARGS.reference_dir / ("reference_" + backend_name + suffix)
        self.assert_ok(run(reference, "--config", native_config, cwd=self.cwd))

        native_files = sorted(p.name for p in ref_package.iterdir())
        self.assertTrue(any(f.startswith("data.") for f in native_files))
        self.assertTrue(any(f.startswith("in.") for f in native_files))
        for filename in native_files:
            with self.subTest(file=filename):
                self.assertTrue((package / filename).is_file())
                # Native info embeds config/output paths; compare its resolved
                # content after mapping those two deliberately different paths.
                expected = (ref_package / filename).read_text()
                actual = (package / filename).read_text()
                if filename.endswith(".info"):
                    def normalize(value, config_path, output_parent):
                        if isinstance(value, dict):
                            return {k: normalize(v, config_path, output_parent) for k, v in value.items()}
                        if isinstance(value, list):
                            return [normalize(v, config_path, output_parent) for v in value]
                        if isinstance(value, str):
                            for spelling in (str(config_path), config_path.as_posix()):
                                value = value.replace(spelling, "CONFIG")
                            for spelling in (str(output_parent), output_parent.as_posix()):
                                value = value.replace(spelling, "OUTPUT")
                        return value
                    self.assertEqual(
                        normalize(json.loads(actual), package / "generator.conf", self.root / "unified"),
                        normalize(json.loads(expected), native_config, ref_parent))
                else:
                    self.assertEqual(actual, expected)

        record = json.loads((package / "siliconelab.json").read_text())
        self.assertEqual(record["system"], system)
        self.assertEqual(record["model"], model)
        self.assertEqual(record["status"], "generated")
        self.assertEqual(record["source_config"], config.as_posix())
        self.assertEqual(record["data_file"], "data." + name)
        self.assertEqual(set(record["generated_files"]), {p.name for p in package.iterdir()})
        self.assertEqual((package / "request.conf").read_bytes(), config.read_bytes())
        source = "network" if system == "elastomer" else system
        upstream = json.loads((ARGS.source_dir / "upstream.json").read_text())
        self.assertEqual(record["upstream"]["commit"], upstream["sources"][source]["integrated_commit"])

    def test_examples_match_native_generators(self):
        for system, model, filename in [
            ("oil", None, "oil.conf"),
            ("elastomer", None, "elastomer.conf"),
            ("coating", "v22", "coating-v22.conf"),
            ("coating", "v35", "coating-v35.conf"),
        ]:
            with self.subTest(system=system, model=model):
                self.parity(system, model, example_entries(filename), filename[:-5])

    def test_repeated_chain_counts_and_copolymer(self):
        self.parity("oil", None, [
            ("chain_count", "8 4"), ("chain_count", "16 4"),
            ("mps_percent", "50"), ("sequence", "alternating"), ("seed", "31415")
        ], "polydisperse")

    def test_elastomer_star_filler_and_conversion(self):
        self.parity("elastomer", None, [
            ("strand_topology", "star"), ("strand_length", "16"),
            ("strand_count", "12"), ("strand_arm_count", "4"),
            ("strand_functionality", "4"), ("crosslinker_length", "16"),
            ("functionality", "4"), ("stoichiometry", "1:1"),
            ("moderator_count", "0"), ("filler_length", "4"),
            ("filler_wt", "10"), ("target_conversion", "80")
        ], "star-filler")

    def test_elastomer_film(self):
        self.parity("elastomer", None, example_entries("elastomer.conf") + [
            ("thickness", "30")
        ], "elastomer-film")

    def test_coating_surface(self):
        self.parity("coating", "v22", example_entries("coating-v22.conf") + [
            ("thickness", "30"), ("surface_relax_steps", "1000"),
            ("surface_production_steps", "2000"), ("surface_sample_every", "100")
        ], "coating-surface")

    def test_invalid_configs_leave_no_package(self):
        cases = [
            ("", "system must be"),
            ("system = unknown\n", "system must be"),
            ("system = coating\n", "Coating requires model"),
            ("system = coating\nmodel = v99\n", "Coating requires model"),
            ("system = oil\nmodel = v22\n", "only supported"),
            ("system = oil\nsystem = coating\n", "Duplicate setting"),
            ("system = oil\noutput_dir = other\n", "Duplicate setting"),
            ("system = oil\noutput = data.test\n", "Use output_dir"),
            ("system = oil\nconfig = other.conf\n", "Reserved setting"),
            ("system = oil\nlength = -1\n", "Error:"),
            ("system = oil\nstrand_length = 8\n", "Unknown option"),
            ("system = oil\nlength 8\n", "must use key = value"),
        ]
        for index, (text, diagnostic) in enumerate(cases):
            with self.subTest(config=text):
                package = self.root / f"invalid-{index}"
                config = self.root / f"invalid-{index}.conf"
                config.write_text(f"output_dir = {package.as_posix()}\n" + text)
                result = run(ARGS.generator, "--config", config, cwd=self.cwd)
                self.assertNotEqual(result.returncode, 0)
                self.assertIn(diagnostic, result.stderr)
                self.assertFalse(package.exists())

    def test_output_required_and_existing_run_preserved(self):
        config = self.root / "oil.conf"
        config.write_text("system = oil\nlength = 4\nchains = 2\n")
        self.assertIn("Missing required output_dir", run(ARGS.generator, "--config", config, cwd=self.cwd).stderr)
        package = self.root / "existing"
        package.mkdir()
        sentinel = package / "keep-me"
        sentinel.write_text("valuable existing data")
        config.write_text(config.read_text() + "output_dir = existing\n")
        result = run(ARGS.generator, "--config", config, cwd=self.cwd)
        self.assertNotEqual(result.returncode, 0)
        self.assertIn("already exists", result.stderr)
        self.assertEqual(sentinel.read_text(), "valuable existing data")

    def test_all_help_routes(self):
        for route in [(), ("oil",), ("elastomer",), ("coating-v22",), ("coating-v35",)]:
            self.assert_ok(run(ARGS.generator, "--help", *route, cwd=self.cwd))

    def test_run_names_must_work_in_lammps_inputs(self):
        config = self.root / "bad-name.conf"
        for name in ("run with spaces", "run$variable", "run;command"):
            with self.subTest(name=name):
                config.write_text(f'system = oil\noutput_dir = "{name}"\nlength = 4\nchains = 2\n')
                result = run(ARGS.generator, "--config", config, cwd=self.cwd)
                self.assertNotEqual(result.returncode, 0)
                self.assertIn("Run directory name", result.stderr)
                self.assertFalse((self.root / name).exists())

    def test_imported_sources_are_unchanged(self):
        upstream = json.loads((ARGS.source_dir / "upstream.json").read_text())
        for source in upstream["sources"].values():
            self.assertEqual(source["license"], "MIT")
            for entry in source["vendored_files"]:
                with self.subTest(path=entry["path"]):
                    actual = hashlib.sha256((ARGS.source_dir / entry["path"]).read_bytes()).hexdigest()
                    self.assertEqual(actual, entry["sha256"])


if __name__ == "__main__":
    unittest.main(argv=[__file__, *UNITTEST_ARGS])
