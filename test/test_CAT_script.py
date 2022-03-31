from script.generate_cannot_link import generate_CAT
import pandas as pd

def test_CAT_script():
    cat_result = pd.read_csv('test/CAT_data/CAT.out', sep='\t')
    cat_result = cat_result[['# contig', 'classification', 'genus', 'species']]
    cat_result = cat_result.values
    _, cannot_link_genus, cannot_link_species = generate_CAT(cat_result)
    assert len(cannot_link_species) == 1
    assert len(cannot_link_genus) == 257

if __name__ == '__main__':
    test_CAT_script()

