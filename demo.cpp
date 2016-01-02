#include "util.h"
#include "detect_keypoints.h"
#include "ddd.h"
#include "fragIO.h"
#include "oldFormatReader.h"

int main(int argc, char **argv) {






  //------------------ Demo for loading old formatted TSDFs ------------------//

  // // Convert TSDF from old format to new format
  // std::string old_fragment_name = "/data/andyz/kinfu/sun3d/mit_76_studyroom/76-1studyroom2/scene50_99";
  // std::string new_fragment_name = "/data/andyz/fragments/sun3d/mit_176_studyroom2_50_99";
  // convert_tsdf_old_to_new(old_fragment_name, new_fragment_name);
  // old_fragment_name = "/data/andyz/kinfu/sun3d/mit_76_studyroom/76-1studyroom2/scene100_149";
  // new_fragment_name = "/data/andyz/fragments/sun3d/mit_176_studyroom2_100_149";
  // convert_tsdf_old_to_new(old_fragment_name, new_fragment_name);
  // old_fragment_name = "/data/andyz/kinfu/train/fire_seq01/scene0_30";
  // new_fragment_name = "/data/andyz/fragments/7scenes/fire_seq01_0_30";
  // convert_tsdf_old_to_new(old_fragment_name, new_fragment_name);
  // old_fragment_name = "/data/andyz/kinfu/train/fire_seq01/scene780_810";
  // new_fragment_name = "/data/andyz/fragments/7scenes/fire_seq01_780_810";
  // convert_tsdf_old_to_new(old_fragment_name, new_fragment_name);

  // Convert point cloud to fragment file (ignore for now)







  //------------------ Demo for aligning 2 fragments ------------------//

  // Get keypoints and keypoint features from first fragment
  std::string frag1_dir = "/data/andyz/fragments/7scenes/fire_seq01_0_30";
  std::cerr << "Loading fragment: " << frag1_dir << std::endl;

  std::vector<std::vector<int>> grid_keypoints1;
  std::vector<std::vector<float>> world_keypoints1;
  get_frag_keypoints(frag1_dir, grid_keypoints1, world_keypoints1);

  std::vector<std::vector<float>> keypoint_features1;
  get_frag_features(frag1_dir, grid_keypoints1, 2048, keypoint_features1);

  // Get keypoints and keypoint features from second fragment
  std::string frag2_dir = "/data/andyz/fragments/7scenes/fire_seq01_780_810";
  std::cerr << "Loading fragment: " << frag2_dir << std::endl;

  std::vector<std::vector<int>> grid_keypoints2;
  std::vector<std::vector<float>> world_keypoints2;
  get_frag_keypoints(frag2_dir, grid_keypoints2, world_keypoints2);

  std::vector<std::vector<float>> keypoint_features2;
  get_frag_features(frag2_dir, grid_keypoints2, 2048, keypoint_features2);

  // Compare keypoint features between two fragments
  std::string frag_cache_dir = "/data/andyz/fragments/cache";
  std::string frag_pair_name = "fire_seq01_0_30_780_810";
  std::vector<std::vector<float>> score_matrix1, score_matrix2;
  compare_frag_features(keypoint_features1, keypoint_features2, frag_cache_dir, frag_pair_name, score_matrix1, score_matrix2);
  save_frag_world_keypoints_cache(world_keypoints1, world_keypoints2, frag_cache_dir, frag_pair_name);

  // // RANSAC parameters
  // float k_match_score_thresh = 0.1f;
  // float ransac_k = 10;
  // float max_ransac_iter = 1000000;
  // float ransac_thresh = 0.04f;

  // // Use RANSAC to compute rigid transform
  // float* Rt = new float[16]; 
  // Rt[12] = 0; Rt[13] = 0; Rt[14] = 0; Rt[15] = 1;
  // std::vector<std::vector<float>> inlier_pairs;
  // ddd_align_feature_cloud(world_keypoints1, keypoint_features1, score_matrix1,
  //                         world_keypoints2, keypoint_features2, score_matrix2,
  //                         k_match_score_thresh, ransac_k, max_ransac_iter, ransac_thresh, Rt, inlier_pairs);

  // // Show computed Rt
  // for (int i = 0; i < 16; i++)
  //   std::cout << Rt[i] << std::endl;
  // for (int i = 0; i < inlier_pairs.size(); i++) {
  //   std::cout << inlier_pairs[i][0] << " " << inlier_pairs[i][1] << " " << inlier_pairs[i][2] << " " << inlier_pairs[i][3] << " " << inlier_pairs[i][4] << " " << inlier_pairs[i][5] << std::endl;
  // }

  return 0;
}