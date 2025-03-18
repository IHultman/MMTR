Steps to create the LFW real data analysis dataset:
  * Download the LFW deepfunneled dataset from https://vis-www.cs.umass.edu/lfw/
  * Run the 'create_lfw_subset_for_mmtr.R' script to filter to people who have between 4 and 50 
    images.
  * Run the 'resize_images.R' script to grayscale and resize the images to 32x32.
  * Run the 'create_multiple_train_test_sets.R' script to compute the first principal component of 
    the attributes for all images included in the analysis and create ten train/test datasets.
