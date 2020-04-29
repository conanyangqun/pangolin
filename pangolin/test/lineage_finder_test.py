import unittest
import dendropy
import os

from pangolin.scripts.lineage_finder import LineageFinder, all_equal, trim_to_common_ancestor, \
    get_basal_lineage
from pangolin.scripts.utils import collapse_nodes

this_dir = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
data_dir = os.path.join(this_dir, "test", 'data', 'lineage_finder')


class LineageTests(unittest.TestCase):
    def test_all_equal(self):
        self.assertTrue(all_equal(["B.1.1", "B.1.1", "B.1.1"]))
        self.assertFalse(all_equal(["B.2.1", "B.1.1"]))

    def test_trim(self):
        self.assertEqual(trim_to_common_ancestor(["B.1.2", "B.1.3", "B.1.4"]), "B.1")
        self.assertEqual(trim_to_common_ancestor(["B.2.6", "B.1.5", "B.3.4"]), "B")

    def test_basal_lineage(self):
        simple = ["B.1.2", "B.1.3", "B.1.4"]
        a_little_tricky = ["A", "B.2"]
        extra_ticky = ["A.1", "A.1.2", "B.3", "C.3", "C.3"]

        self.assertEqual(get_basal_lineage(simple), "B.1")
        self.assertEqual(get_basal_lineage(a_little_tricky), "A")
        self.assertEqual(get_basal_lineage(extra_ticky), "A.1")

    def test_sibling(self):
        tree = dendropy.Tree.get_from_string("(A|A,(B|B,((C|B.1,test),D|B.1)77));", "newick")
        finder = LineageFinder(tree, "test", 1, "|")

        self.assertEqual(finder.get_lineage(), ["B.1", "77"])

    def test_nested_no_sibling(self):
        tree = dendropy.Tree.get_from_string("(A|A,(B|B,((C|B.1,c1|B.1),(D|B.1,d1|B.1),test)77));", "newick")

        finder = LineageFinder(tree, "test", 1, "|")

        self.assertEqual(finder.get_lineage(), ["B.1", "77"])

    def test_between_clades_fall_back(self):
        tree = dendropy.Tree.get_from_string("(A|A,((B|B,(C|B.1,D|B.1)),test))77;", "newick")

        finder = LineageFinder(tree, "test", 1, "|")

        self.assertEqual(finder.get_lineage(), ["A", "77"])

    def test_polytomy_no_tip_Siblings(self):
        tree = dendropy.Tree.get_from_string("(A|A,((B|B,b1|B),test,(C|B.1,D|B.1))77,a1|A);", "newick")
        finder = LineageFinder(tree, "test", 1, "|")

        self.assertEqual(finder.get_lineage(), ["B", "77"])


if __name__ == '__main__':
    unittest.main()
