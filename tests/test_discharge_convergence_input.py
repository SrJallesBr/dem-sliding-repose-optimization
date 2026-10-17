import re
import unittest
from pathlib import Path


class DischargeDurationTests(unittest.TestCase):
    def test_discharge_uses_15_second_fixed_duration_and_writes_final_snapshot(self):
        deck = (Path(__file__).resolve().parents[1] / "Ang Deslizamento" / "in.deslizamento").read_text(encoding="utf-8")
        self.assertRegex(deck, r"(?m)^variable\s+simulationtime\s+equal\s+15(?:\.0)?\s*(?:#.*)?$")
        self.assertIn("variable steps_finais equal round((v_simulationtime-v_tempo_remocao_base)/v_dt)", deck)
        self.assertRegex(deck, r"(?m)^run\s+\$\{steps_finais\}.*$\n\nwrite_dump all custom/vtk post/deslizamento_final\.vtk")


if __name__ == "__main__":
    unittest.main()
