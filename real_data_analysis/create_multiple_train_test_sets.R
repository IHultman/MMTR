library(dplyr)
library(R.matlab)
library(stringr)


data_dir = "./resized_images_for_mmtr/";
lfw_attributes_filename = "./lfw_attributes.csv";
mmtr_save_dir = "./mmtr_lfw_multiple_train_test_sets/";
kruskal_save_dir = "./kruskal_lfw_multiple_train_test_sets/";

if (!dir.exists(mmtr_save_dir) ) {
  dir.create(mmtr_save_dir, recursive=TRUE);

  stopifnot(dir.exists(mmtr_save_dir) );
}

if (!dir.exists(kruskal_save_dir) ) {
  dir.create(kruskal_save_dir, recursive=TRUE);

  stopifnot(dir.exists(kruskal_save_dir) );
}

expected_img_dim = c(32, 32);

ppl_info = merge(
  bind_rows(
    lapply(
      list.files(data_dir),
      \(person_px) {
        person_dir_px = paste0(data_dir, "/", person_px);
        img_filenames_px = list.files(person_dir_px);
        file_id_re = sprintf("(?<=^%s_)[0-9]+(?=\\.RData$)", person_px);
        data.frame(
          person=str_replace_all(person_px, "_", " "),
          imagenum=as.numeric(str_extract(img_filenames_px, file_id_re) ),
          img_filename=paste0(person_dir_px, "/", img_filenames_px) )}) ),
  read.csv(lfw_attributes_filename, header=TRUE),
  by=c("person", "imagenum") );

names_to_keep = (
  ppl_info %>%
  group_by(person) %>%
  summarize(n=n() ) %>%
  filter(n > 3) %>%
  select(person) );

ppl_info = merge(ppl_info, names_to_keep, by="person");
order_ixs = order(ppl_info$person, ppl_info$imagenum);
ppl_info = ppl_info[order_ixs,];
rownames(ppl_info) = NULL;

attributes_pca = prcomp(
  as.matrix(select(ppl_info, !c(person, imagenum, img_filename) )),
  center=TRUE,
  scale.=TRUE);

ppl_info$PC1 = attributes_pca$x[,1];

n_images = nrow(ppl_info);
all_images_arr = array(NA, dim=c(expected_img_dim, n_images) );

for (ix in 1:n_images) {
  img_ix = readRDS(ppl_info$img_filename[ix]);

  stopifnot(all(dim(img_ix) == expected_img_dim) );

  all_images_arr[,,ix] = img_ix;
}

set.seed(13);

n_tst_imgs_per_grp = 4;
min_images_in_test_groups = 12;
n_train_test_sets = 10;

test_groups_ix_ranges = (
  ppl_info %>%
  mutate(rx=1:n_images) %>%
  group_by(person) %>%
  summarize(n=n(), rx_min=min(rx), rx_max=max(rx) ) %>%
  filter(n >= min_images_in_test_groups) %>%
  select(rx_min, rx_max) %>%
  as.matrix() );

for (sx in 1:n_train_test_sets) {
  test_ixs = sort(
    as.numeric(
      apply(
        test_groups_ix_ranges,
        MARGIN=1,
        \(row_rx) sample(row_rx[1]:row_rx[2], n_tst_imgs_per_grp) )));

  train_ixs = setdiff(1:n_images, test_ixs);

  X_train = all_images_arr[,,train_ixs];
  X_test = all_images_arr[,,test_ixs];

  y_train = select(ppl_info[train_ixs,], !c(person, imagenum, img_filename) );
  y_test = select(ppl_info[test_ixs,], !c(person, imagenum, img_filename) );

  group_labels_train = as.factor(ppl_info$person[train_ixs]);
  group_labels_test = as.factor(ppl_info$person[test_ixs]);

  stopifnot(dim(X_train)[3] == nrow(y_train) );
  stopifnot(length(group_labels_train) == nrow(y_train) );
  stopifnot(dim(X_test)[3] == nrow(y_test) );
  stopifnot(length(group_labels_test) == nrow(y_test) );

  save_filename = paste0("/lfw_train_test_set_", sx);

  saveRDS(
    list(
      X_train=X_train,
      X_test=X_test,
      y_train=y_train,
      y_test=y_test,
      group_labels_train=group_labels_train,
      group_labels_test=group_labels_test),
    paste0(mmtr_save_dir, save_filename, ".RData") );

  writeMat(
    paste0(kruskal_save_dir, save_filename, ".mat"),
    X_train=X_train,
    X_test=X_test,
    y_train=y_train,
    y_test=y_test,
    group_labels_train=as.character(group_labels_train),
    group_labels_test=as.character(group_labels_test) );
}
