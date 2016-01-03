// #include "pc2tsdf/pc2tsdf.h"
// #include "detect_keypoints.h"
// #include "ddd.h"
// #include "fragmentMatcher/fragmentMatcher.h"
#include <vector>
#include <numeric>
#include <algorithm>
#include <iostream>
#include <fstream>
#include <stdio.h>
#include <stdlib.h>
#include <string>
#include <cmath>
#include "cluster_ransacK.cpp"
#include <functional>

bool fileExists(const std::string& filename) {
  std::ifstream file(filename.c_str());
  return (!file.fail());
}

bool sort_arr_desc_compare(int a, int b, float* data) {
  return data[a] > data[b];
}

std::vector<std::vector<float>> ddd_align_feature_cloud(const std::vector<std::vector<float>> &world_keypoints1, const std::vector<std::vector<float>> &score_matrix1,
                             const std::vector<std::vector<float>> &world_keypoints2, const std::vector<std::vector<float>> &score_matrix2,
                             float voxelSize, float k_match_score_thresh, float ransac_k, float max_ransac_iter, float ransac_thresh, float* Rt) {

  // For each keypoint from first set, find indices of all keypoints
  // in second set with score > k_match_score_thresh
  std::vector<std::vector<int>> match_rank1;
  for (int i = 0; i < world_keypoints1.size(); i++) {
    // Sort score vector in descending fashion
    std::vector<float> tmp_score_vect = score_matrix1[i];
    float* tmp_score_vect_arr = &tmp_score_vect[0];
    int* tmp_score_idx = new int[tmp_score_vect.size()];
    std::iota(tmp_score_idx, tmp_score_idx + tmp_score_vect.size(), 0);
    std::sort(tmp_score_idx, tmp_score_idx + tmp_score_vect.size(), std::bind(sort_arr_desc_compare, std::placeholders::_1, std::placeholders::_2, tmp_score_vect_arr));
    std::vector<int> tmp_score_rank;
    for (int j = 0; j < world_keypoints2.size(); j++)
      if (tmp_score_vect_arr[tmp_score_idx[j]] > k_match_score_thresh)
        tmp_score_rank.push_back(tmp_score_idx[j]);
    // std::cout << tmp_score_rank.size() << std::endl;
    match_rank1.push_back(tmp_score_rank);
  }

  // For each keypoint from second set, find indices of all keypoints
  // in first set with score > k_match_score_thresh
  std::vector<std::vector<int>> match_rank2;
  for (int i = 0; i < world_keypoints2.size(); i++) {
    // Sort score vector in descending fashion
    std::vector<float> tmp_score_vect = score_matrix2[i];
    float* tmp_score_vect_arr = &tmp_score_vect[0];
    int* tmp_score_idx = new int[tmp_score_vect.size()];
    std::iota(tmp_score_idx, tmp_score_idx + tmp_score_vect.size(), 0);
    std::sort(tmp_score_idx, tmp_score_idx + tmp_score_vect.size(), std::bind(sort_arr_desc_compare, std::placeholders::_1, std::placeholders::_2, tmp_score_vect_arr));
    std::vector<int> tmp_score_rank;
    for (int j = 0; j < world_keypoints1.size(); j++)
      if (tmp_score_vect_arr[tmp_score_idx[j]] > k_match_score_thresh)
        tmp_score_rank.push_back(tmp_score_idx[j]);
    // std::cout << tmp_score_rank.size() << std::endl;
    match_rank2.push_back(tmp_score_rank);
  }

  // Finalize match matrix (indices) unofficial reflexive property
  // A pair of points (with feature vectors f1 and f2) match iff
  // ddd(f1,f2) > threshold && ddd(f2,f1) > threshold
  std::vector<std::vector<int>> match_idx;
  for (int i = 0; i < world_keypoints1.size(); i++) {
    std::vector<int> tmp_matches;
    for (int j = 0; j < match_rank1[i].size(); j++) {
      int tmp_match_idx = match_rank1[i][j];
      if (std::find(match_rank2[tmp_match_idx].begin(), match_rank2[tmp_match_idx].end(), i) != match_rank2[tmp_match_idx].end())
        tmp_matches.push_back(tmp_match_idx);
    }
    match_idx.push_back(tmp_matches);
  }

  // DEBUG
  for (int i = 0; i < world_keypoints1.size(); i++) {
    std::cout << i << " | ";
    for (int j = 0; j < match_idx[i].size(); j++)
      std::cout << match_idx[i][j] << " ";
    std::cout << std::endl;
  }

  // Compute Rt transform from second to first point cloud (k-ransac)
  std::vector<std::vector<float>> inliers = ransacfitRt(world_keypoints1, world_keypoints2, match_idx, ransac_k, max_ransac_iter, ransac_thresh, Rt, true);
  return inliers;
}

void loadKeypoints(std::string &filename, std::vector<std::vector<float>> &grid) {

  if (!fileExists(filename)) {
    std::cout << "Could not read file: " << filename << std::endl;
    exit(1);
  }

  FILE *fp = fopen(filename.c_str(), "r");
  int num_keypoints = 0;
  int iret = fscanf(fp, "# number of keypoints: %d", &num_keypoints);
  for (int i = 0; i < num_keypoints; i++) {
    std::vector<float> tmp_world_keypoint;
    for (int j = 0; j < 3; j++) {
      float tmp_coord;
      iret = fscanf(fp, "%f", &tmp_coord);
      tmp_world_keypoint.push_back(tmp_coord);
    }
    grid.push_back(tmp_world_keypoint);
  }
  fclose(fp);
}

void loadScores(std::string &filename, std::vector<std::vector<float>> &grid) {
  
  if (!fileExists(filename)) {
    std::cout << "Could not read file: " << filename << std::endl;
    exit(1);
  }

  FILE *fp = fopen(filename.c_str(), "r");
  int num_rows = 0;
  int num_cols = 0;
  int iret = fscanf(fp, "# number of keypoint comparisons: %d x %d", &num_rows, &num_cols);
  for (int i = 0; i < num_rows; i++) {
    std::vector<float> tmp_scores;
    for (int j = 0; j < num_cols; j++) {
      float tmp_coord;
      iret = fscanf(fp, "%f", &tmp_coord);
      tmp_scores.push_back(tmp_coord);
    }
    grid.push_back(tmp_scores);
  }
  fclose(fp);
}


int main(int argc, char **argv) {

  // ./main  frag_pair_name 

  std::string frag_cache_dir = "/n/fs/sun3d/DDD/cache/";
  std::string frag_pair_name = argv[1];

  // Load first fragment's keypoints
  std::vector<std::vector<float>> keypoints1;
  std::string keypoints1_filename = frag_cache_dir + "/" + frag_pair_name + "_keypoints1.txt";
  loadKeypoints(keypoints1_filename, keypoints1);

  // Load second fragment's keypoints
  std::vector<std::vector<float>> keypoints2;
  std::string keypoints2_filename = frag_cache_dir + "/" + frag_pair_name + "_keypoints2.txt";
  loadKeypoints(keypoints2_filename, keypoints2);

  // Load first score matrix
  std::vector<std::vector<float>> score_matrix1;
  std::string score_matrix1_filename = frag_cache_dir + "/" + frag_pair_name + "_scores1.txt";
  loadScores(score_matrix1_filename, score_matrix1);

  // Load second score matrix
  std::vector<std::vector<float>> score_matrix2;
  std::string score_matrix2_filename = frag_cache_dir + "/" + frag_pair_name + "_scores2.txt";
  loadScores(score_matrix2_filename, score_matrix2);

  // Contains rigid transform matrix
  float* Rt = new float[16];
  Rt[12] = 0; Rt[13] = 0; Rt[14] = 0; Rt[15] = 1;

  // RANSAC parameters
  const float k_match_score_thresh = 0.1f;
  const float ransac_k = 10; // RANSAC over top-k > k_match_score_thresh
  const float max_ransac_iter = 2000000;
  const float ransac_inlier_thresh = 0.04f;
  const float voxelSize = 0.01f;

  std::vector<std::vector<float>> inliers = ddd_align_feature_cloud(keypoints1, score_matrix1, keypoints2, score_matrix2, voxelSize, k_match_score_thresh, ransac_k, max_ransac_iter, ransac_inlier_thresh, Rt);

  // Save to file in results
  std::string rt_filename = frag_cache_dir + "/" + frag_pair_name + "_rt.txt";
  FILE *file = fopen(rt_filename.c_str(), "wb");
  for (int i = 0; i < 4; i++)
    fprintf(file, "%.17g %.17g %.17g %.17g\n", Rt[4 * i + 0], Rt[4 * i + 1], Rt[4 * i + 2], Rt[4 * i + 3]);

  fprintf(file, "num_inliers: %d\n",inliers.size());
  fclose(file);


  return 0;
}