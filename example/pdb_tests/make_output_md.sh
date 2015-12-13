#!/bin/bash
# Generates notebook in markdown
## pip install jupyter notedown

# markdown to markdown with outputs
notedown _input_example-mapping_and_distances.md --run --render --timeout=600 |\
    perl -nle 's/<.?div.*$/```/ ; print' > MAPPING_AND_DISTANCES.md

# markdown to jupyter
#notedown _input_example-mapping_and_distances.md > _jupyter_example-mapping_and_distance.ipynb

# jupyter to markdown with outputs
#notedown _jupyter_example-mapping_and_distance.ipynb --render |\
#    perl -nle 's/<.?div.*$/```/ ; print' > MAPPING_AND_DISTANCES.md
