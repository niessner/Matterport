# Semantic annotation metadata

The object instance-level semantic annotations are specified by annotators as freeform text, which we then post-process and correspond to a variety of category taxonomies.  We provide the following files:

- [category_mapping.tsv](category_mapping.tsv) : A mapping of all raw object category text labels to category taxonomies such as WordNet, NYUv2, ShapeNet, and our own 40 category set (mp40cat)
- [mpcat40.tsv](mpcat40.tsv) : Metadata for each of our 40 categories, including visualization color ("hex" column), WordNet synset key, corresponding NYUv2 40 label, and labels which are mapped to this category
