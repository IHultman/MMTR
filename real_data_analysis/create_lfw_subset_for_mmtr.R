library(dplyr)


lfw_data_dir = "./lfw-deepfunneled/";
mmtr_subset_save_dir = "./lfw-deepfunneled_for_mmtr/";

if (!dir.exists(mmtr_subset_save_dir) ) {
  dir.create(mmtr_subset_save_dir, recursive=TRUE);

  stopifnot(dir.exists(mmtr_subset_save_dir) );
}

# Make table matching each name with their respective image files.
ppl_file_info = bind_rows(
  lapply(
    list.files(lfw_data_dir),
    \(person_px) {
      person_dir_px = paste0(lfw_data_dir, "/", person_px);
      data.frame(
        name=person_px,
        img_filename=list.files(person_dir_px) )}) );

# Find names who have a number of images within the specified range.
names_to_keep = (
  ppl_file_info %>%
  group_by(name) %>%
  summarize(n=n() ) %>%
  filter(n > 3 & n < 51) %>%
  select(name) );

for (px in 1:nrow(names_to_keep) ) {
  person_px = names_to_keep$name[px];
  from_dir_px = paste0(lfw_data_dir, "/", person_px);
  to_dir_px = paste0(mmtr_subset_save_dir, "/", person_px);
  R.utils::copyDirectory(from_dir_px, to_dir_px);
}
