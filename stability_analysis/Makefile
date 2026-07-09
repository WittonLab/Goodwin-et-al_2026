.SILENT:
preprocessing:
	echo "Concatenating coordinates..."
	python3 utils/session_preprocessing.py $(DATASET)
	echo "Extracting laps from coordinates..."
	python3 utils/pull_laps.py $(DATASET)
	echo "Concatenating sessions..."
	python3 utils/concat_sessions.py $(DATASET)
	echo "Finished."
