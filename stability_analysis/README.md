# Setup and Use
To start, clone this repo locally. Any identically formatted data can still be run through this pipeline.


`git clone https://github.com/WittonLab/Goodwin-et-al_2026.git`

Create a python virtual environment like so:

`python3 -m venv venv`

Install the packages from the requirements file:

`python3 -m pip install -r requirements.txt`

To activate your virtual environment, from the root of the repo, run:

`source venv/bin/activate`

The script `pull_laps.py` generates the lap starts for all mice and experiments. `session_preprocessing.py` pulls out the behaviours and aligns them to the calcium traces. `concat_sessions.py` concatenates the pre-switch recording to the post-switch recording.

`main.py` calls the relevant manifold and RQA functions.

There is a Makefile included which runs the preprocessing steps. Run the Makefile before calling the main script:

`make preprocessing DATASET=XXX`,

where XXX is the path to the dataset containing the calcium traces. 

After running the Makefile, call the main script to start analysis:

`python3 main.py XXX`,

This may take a few minutes to complete. Calling `main.py` will generate the manifold and correlation plots.

`ach_behav_analysis` analyses another cohort of mice from which acetylcholine signals were recorded. This script runs manifold analysis on the behavioural recordings of these mice and correlates it with the levels of acetylcholine.

