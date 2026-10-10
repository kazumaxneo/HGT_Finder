import json
import tempfile
import unittest
from pathlib import Path
from unittest.mock import patch

from gtdb_hgt.cli import parser
from gtdb_hgt.search import SENSITIVITY_MODES, search


class SearchModeTests(unittest.TestCase):
    def test_omitted_mode_uses_very_sensitive(self):
        for command, inputs in (
            ("search", ["--query-faa", "q.faa"]),
            ("run", ["--genbank", "q.gbk", "--taxonomy", "taxonomy.tsv", "--query-genome", "GCF_000000001.1"]),
        ):
            with self.subTest(command=command):
                args = parser().parse_args([command, *inputs, "--db", "db.dmnd", "--out", "output"])
                self.assertEqual(args.sensitivity, "very-sensitive")
        with tempfile.TemporaryDirectory() as d:
            output = Path(d) / "hits.tsv"
            with patch("gtdb_hgt.search.execute") as execute:
                manifest = search(Path("q.faa"), Path("db.dmnd"), output, "diamond", 4)
            self.assertIn("--very-sensitive", execute.call_args.args[0])
            self.assertEqual(manifest["sensitivity"], "very-sensitive")

    def test_cli_modes_and_command_provenance(self):
        with tempfile.TemporaryDirectory() as d:
            root = Path(d)
            for mode in SENSITIVITY_MODES:
                with self.subTest(mode=mode):
                    output = root / (mode + ".tsv.gz")
                    args = parser().parse_args(["search", "--query-faa", "q.faa", "--db", "db.dmnd", "--out", str(output), "--sensitivity", mode])
                    self.assertEqual(args.sensitivity, mode)
                    with patch("gtdb_hgt.search.execute") as execute:
                        search(root / "q.faa", root / "db.dmnd", output, "diamond", 4, sensitivity=mode)
                    cmd = execute.call_args.args[0]
                    selected_flags = ["--" + m for m in SENSITIVITY_MODES if "--" + m in cmd]
                    self.assertEqual(selected_flags, [] if mode == "default" else ["--" + mode])
                    self.assertEqual(cmd[cmd.index("--max-target-seqs") + 1], "0")
                    self.assertIn("pident", cmd)
                    self.assertIn("--compress", cmd)
                    manifest = json.loads(output.with_suffix(".gz.json").read_text())
                    self.assertEqual(manifest["sensitivity"], mode)
                    self.assertEqual(manifest["command"], cmd)

    def test_invalid_mode_fails_before_writing(self):
        with tempfile.TemporaryDirectory() as d:
            output = Path(d) / "new" / "hits.tsv"
            with patch("gtdb_hgt.search.execute") as execute:
                with self.assertRaisesRegex(ValueError, "Unknown DIAMOND sensitivity"):
                    search(Path("q"), Path("db"), output, "diamond", 4, sensitivity="ultra-fast")
                execute.assert_not_called()
            self.assertFalse(output.parent.exists())
