import re
import unittest
from pathlib import Path


PROJECT_ROOT = Path(__file__).resolve().parents[1]
MAX_OBSERVED_RELATIVE_TRAVEL_M = 0.006849
SAFETY_FACTOR = 2.0
TEST_DECKS = (
    PROJECT_ROOT / "Ang Deslizamento" / "in.deslizamento",
    PROJECT_ROOT / "Ang Repouso" / "in.repouso",
)


class NeighborSkinSafetyTests(unittest.TestCase):
    def test_neighbor_skin_exceeds_twice_maximum_observed_relative_travel(self):
        for deck_path in TEST_DECKS:
            with self.subTest(deck=deck_path.parent.name):
                deck = deck_path.read_text(encoding="utf-8")
                match = re.search(
                    r"(?m)^\s*variable\s+neighborskin\s+equal\s+([0-9.eE+-]+)",
                    deck,
                )
                self.assertIsNotNone(match, f"neighborskin variable missing in {deck_path}")
                skin_m = float(match.group(1))
                required_m = MAX_OBSERVED_RELATIVE_TRAVEL_M * SAFETY_FACTOR
                self.assertGreaterEqual(
                    skin_m,
                    required_m,
                    f"{deck_path} skin {skin_m:g} m is below 2x observed travel {required_m:g} m",
                )


if __name__ == "__main__":
    unittest.main()
