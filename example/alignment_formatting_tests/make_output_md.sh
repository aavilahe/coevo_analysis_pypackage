#!/bin/bash
# Generates notebook in markdown
## pip install jupyter notedown

# markdown to markdown with outputs
notedown _input_example-sequence_formatting.md --run --render |\
    perl -nle 's/<.?div.*$/```/ ; print' > SEQUENCE_FORMATTING.md

# markdown to jupyter
#notedown _input_example-sequence_formatting.md > _jupyter_example-mapping_and_distance.ipynb

# jupyter to markdown with outputs
#notedown _jupyter_example-mapping_and_distance.ipynb --render |\
#    perl -nle 's/<.?div.*$/```/ ; print' > SEQUENCE_FORMATTING.md
