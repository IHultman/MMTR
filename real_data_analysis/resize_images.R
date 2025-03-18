library(dplyr)
library(imager)
library(stringr)


lfw_data_dir = "./lfw-deepfunneled_for_mmtr/";
resized_images_save_dir = "./resized_images_for_mmtr/";

if (!dir.exists(resized_images_save_dir) ) {
  dir.create(resized_images_save_dir, recursive=TRUE);

  stopifnot(dir.exists(resized_images_save_dir) );
}

ppl_names = list.files(lfw_data_dir);

resize_dims = c(32, 32);

for (px in 1:length(ppl_names) ) {
  person_px = ppl_names[px];
  from_dir_px = paste0(lfw_data_dir, "/", person_px);
  to_dir_px = paste0(resized_images_save_dir, "/", person_px);

  if (!dir.exists(to_dir_px) ) {
    dir.create(to_dir_px);

    stopifnot(dir.exists(to_dir_px) );
  }

  img_filenames_px = list.files(from_dir_px);
  n_images_px = length(img_filenames_px);

  for (ix in 1:n_images_px) {
    paste0(from_dir_px, "/", img_filenames_px[ix]) %>%
    load.image() %>%
    grayscale() %>%
    resize(size_x=resize_dims[1], size_y=resize_dims[2]) %>%
    as.matrix() %>%
    t() %>%
    saveRDS(
      sprintf(
        "%s/%s.RData",
        to_dir_px,
        tools::file_path_sans_ext(img_filenames_px[ix]) ));
  }
}
