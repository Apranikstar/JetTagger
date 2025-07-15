import uproot
import matplotlib.pyplot as plt
import numpy as np

# Open the ROOT file
file = uproot.open("test_output.root")
events = file["Events;1"]

# Extract arrays
scores = events["score_recojet_isT"].array()
labels = events["_label_"].array()

# Separate scores by true label
signal_scores = scores[labels == 1]
background_scores = scores[labels == 0]

# Plot histograms
plt.figure(figsize=(8,6))
plt.hist(signal_scores, bins=50, alpha=0.6, label='Signal (true top jets)')
plt.hist(background_scores, bins=50, alpha=0.6, label='Background')
plt.xlabel("Predicted score for isT")
plt.ylabel("Number of jets")
plt.title("Distribution of Predicted Scores")
plt.legend()
plt.grid(True)
plt.tight_layout()
plt.savefig("score_distribution.png")  # Save to file
plt.show()
