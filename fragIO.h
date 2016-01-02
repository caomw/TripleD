

// TSDF Format:
// 3 ints (Xdim,Ydim,Zdim = volume size), followed by float array A of TSDF
// A[z*Ydim*Xdim + y*Xdim + x] = tsdf(x,y,z)
// Note: x,y,z starts at 0,0,0

////////////////////////////////////////////////////////////////////////////////

std::vector<int> load_frag_dim(const std::string &frag_dir) {

  std::string filename = frag_dir + "/volume.tsdf";
  std::ifstream in(filename, std::ios::binary | std::ios::in);

  // Read TSDF volume dimensions
  std::vector<int> dim;
  for (int i = 0; i < 3; i++) {
    int td;
    in.read((char*)&td, sizeof(int));
    dim.push_back(td);
  }
  in.close();
  return dim;
}

////////////////////////////////////////////////////////////////////////////////

void load_frag_tsdf(const std::string &frag_dir, float* tsdf, int x_dim, int y_dim, int z_dim) {

  std::string filename = frag_dir + "/volume.tsdf";
  std::ifstream in(filename, std::ios::binary | std::ios::in);

  // Check TSDF volume dimensions
  std::vector<int> dim;
  for (int i = 0; i < 3; i++) {
    int td;
    in.read((char*)&td, sizeof(int));
    dim.push_back(td);
  }
  if ((x_dim != dim[0]) || (y_dim != dim[1]) || (z_dim != dim[2])) {
    std::cerr << "Error loading fragment " + frag_dir + ": inconsistent dimensions!" << std::endl;
    return;
  }

  // Load TSDF values into volume
  for (int i = 0; i < x_dim*y_dim*z_dim; i++)
    in.read((char*)&tsdf[i], sizeof(float));
  in.close();
}

////////////////////////////////////////////////////////////////////////////////

void save_frag_tsdf(const std::string &frag_dir, float* tsdf, int x_dim, int y_dim, int z_dim) {

  std::string filename = frag_dir + "/volume.tsdf";
  std::ofstream out(filename, std::ios::binary | std::ios::out);

  // Write TSDF volume dimensions
  out.write((char*)&x_dim, sizeof(int));
  out.write((char*)&y_dim, sizeof(int));
  out.write((char*)&z_dim, sizeof(int));

  // Save TSDF values to file
  out.write((char*)tsdf, x_dim * y_dim * z_dim * sizeof(float));
  out.close();
}

////////////////////////////////////////////////////////////////////////////////

void load_frag_ext(const std::string &frag_dir, float* grid2cam, float* cam2world) {

  // Read matrix for converting voxel grid coordinates to camera coordinates
  std::string filename = frag_dir + "/grid2cam.txt";
  FILE *fp = fopen(filename.c_str(), "r");
  for (int i = 0; i < 16; i++)
    int iret = fscanf(fp, "%f", &grid2cam[i]);
  fclose(fp);

  // Read matrix for converting camera coordinates to world coordinates
  filename = frag_dir + "/cam2world.txt";
  fp = fopen(filename.c_str(), "r");
  for (int i = 0; i < 16; i++)
    int iret = fscanf(fp, "%f", &cam2world[i]);
  fclose(fp);
}

////////////////////////////////////////////////////////////////////////////////

void tsdf2ply(const std::string &filename, float* tsdf, float threshold, int x_dim, int y_dim, int z_dim, float* transform) {


  // Count total number of points in point cloud
  int num_points = 0;
  for (int i = 0; i < x_dim * y_dim * z_dim; i++)
    if (std::abs(tsdf[i]) < threshold)
      num_points++;

  // Create header for ply file
  FILE *fp = fopen(filename.c_str(), "w");
  fprintf(fp, "ply\n");
  fprintf(fp, "format binary_little_endian 1.0\n");
  fprintf(fp, "element vertex %d\n", num_points);
  fprintf(fp, "property float x\n");
  fprintf(fp, "property float y\n");
  fprintf(fp, "property float z\n");
  // fprintf(fp, "property uchar red\n");
  // fprintf(fp, "property uchar green\n");
  // fprintf(fp, "property uchar blue\n");
  fprintf(fp, "end_header\n");


  // for (int i = 0; i < 16; i++)
  //   std::cout << transform[i] << std::endl;

  // Create point cloud content for ply file
  for (int z = 0; z < z_dim; z++) {
    for (int y = 0; y < y_dim; y++) {
      for (int x = 0; x < x_dim; x++) {
        if (std::abs(tsdf[z * y_dim * x_dim + y * x_dim + x]) < threshold) {

          // grid to world coords (7scene)
          // float sx = ((float)x + 1) * 0.01 - 512 * 0.01 / 2;
          // float sy = ((float)y + 1) * 0.01 - 512 * 0.01 / 2;
          // float sz = ((float)z + 1) * 0.01 - 0.5;

          // (synth)
          // float sx = (float)x * 0.01;
          // float sy = (float)y * 0.01;
          // float sz = (float)z * 0.01;


          float sx = (float)x;
          float sy = (float)y;
          float sz = (float)z;


          float fx = transform[0] * sx + transform[1] * sy + transform[2] * sz + transform[3];
          float fy = transform[4] * sx + transform[5] * sy + transform[6] * sz + transform[7];
          float fz = transform[8] * sx + transform[9] * sy + transform[10] * sz + transform[11];
          fwrite(&fx, sizeof(float), 1, fp);
          fwrite(&fy, sizeof(float), 1, fp);
          fwrite(&fz, sizeof(float), 1, fp);

          // if (use_ext) {
          //   float sx = (float)x;
          //   float sy = (float)y;
          //   float sz = (float)z;
          //   float fx = transform[0] * sx + transform[1] * sy + transform[2] * sz + transform[3];
          //   float fy = transform[4] * sx + transform[5] * sy + transform[6] * sz + transform[7];
          //   float fz = transform[8] * sx + transform[9] * sy + transform[10] * sz + transform[11];
          //   fx = fx + transform[3];
          //   fy = fy + transform[7];
          //   fz = fz + transform[11];
          //   fwrite(&fx, sizeof(float), 1, fp);
          //   fwrite(&fy, sizeof(float), 1, fp);
          //   fwrite(&fz, sizeof(float), 1, fp);

          //   uchar r = (uchar) pc_color.x;
          //   uchar g = (uchar) pc_color.y;
          //   uchar b = (uchar) pc_color.z;
          //   fwrite(&r, sizeof(uchar), 1, fp);
          //   fwrite(&g, sizeof(uchar), 1, fp);
          //   fwrite(&b, sizeof(uchar), 1, fp);
          // } else {
          //   float sx = (float)x;
          //   float sy = (float)y;
          //   float sz = (float)z;
          //   fwrite(&sx, sizeof(float), 1, fp);
          //   fwrite(&sy, sizeof(float), 1, fp);
          //   fwrite(&sz, sizeof(float), 1, fp);

          //   uchar r = (uchar) pc_color.x;
          //   uchar g = (uchar) pc_color.y;
          //   uchar b = (uchar) pc_color.z;
          //   fwrite(&r, sizeof(uchar), 1, fp);
          //   fwrite(&g, sizeof(uchar), 1, fp);
          //   fwrite(&b, sizeof(uchar), 1, fp);
          // }


          // transform
          // float fx = ext_mat[0] * sx + ext_mat[1] * sy + ext_mat[2] * sz;
          // float fy = ext_mat[4] * sx + ext_mat[5] * sy + ext_mat[6] * sz;
          // float fz = ext_mat[8] * sx + ext_mat[9] * sy + ext_mat[10] * sz;
          // fx = fx + ext_mat[3];
          // fy = fy + ext_mat[7];
          // fz = fz + ext_mat[11];
          // fwrite(&fx, sizeof(float), 1, fp);
          // fwrite(&fy, sizeof(float), 1, fp);
          // fwrite(&fz, sizeof(float), 1, fp);
        }
      }
    }
  }
  fclose(fp);

}

////////////////////////////////////////////////////////////////////////////////

// Load frag (in world coordinates), apply transform, and save to point cloud
void frag2ply(const std::string &frag_dir, const std::string &filename, float surface_tsdf_threshold, float* transform, std::vector<float> &color) {

  // Load TSDF volume
  std::vector<int> tsdf_dim = load_frag_dim(frag_dir);
  float *tsdf = new float[tsdf_dim[0] * tsdf_dim[1] * tsdf_dim[2]];
  load_frag_tsdf(frag_dir, tsdf, tsdf_dim[0], tsdf_dim[1], tsdf_dim[2]);

  // Count total number of points in point cloud
  int num_points = 0;
  for (int i = 0; i < tsdf_dim[0] * tsdf_dim[1] * tsdf_dim[2]; i++)
    if (std::abs(tsdf[i]) < surface_tsdf_threshold)
      num_points++;

  // Create header for ply file
  FILE *fp = fopen(filename.c_str(), "w");
  fprintf(fp, "ply\n");
  fprintf(fp, "format binary_little_endian 1.0\n");
  fprintf(fp, "element vertex %d\n", num_points);
  fprintf(fp, "property float x\n");
  fprintf(fp, "property float y\n");
  fprintf(fp, "property float z\n");
  fprintf(fp, "property uchar red\n");
  fprintf(fp, "property uchar green\n");
  fprintf(fp, "property uchar blue\n");
  fprintf(fp, "end_header\n");

  // Load grid to world matrix
  float *grid2cam = new float[16];
  float *cam2world = new float[16];
  float *grid2world = new float[16];
  load_frag_ext(frag_dir, grid2cam, cam2world);
  multiply_matrix(cam2world, grid2cam, grid2world);

  // Create point cloud content for ply file
  for (int z = 0; z < tsdf_dim[2]; z++) {
    for (int y = 0; y < tsdf_dim[1]; y++) {
      for (int x = 0; x < tsdf_dim[0]; x++) {
        if (std::abs(tsdf[z * tsdf_dim[1] * tsdf_dim[0] + y * tsdf_dim[0] + x]) < surface_tsdf_threshold) {

          // Convert grid to world coordinates
          float sx = (float)x;
          float sy = (float)y;
          float sz = (float)z;
          float tx = grid2world[0] * sx + grid2world[1] * sy + grid2world[2] * sz + grid2world[3];
          float ty = grid2world[4] * sx + grid2world[5] * sy + grid2world[6] * sz + grid2world[7];
          float tz = grid2world[8] * sx + grid2world[9] * sy + grid2world[10] * sz + grid2world[11];

          // Apply Rt transform to points
          float fx = transform[0] * tx + transform[1] * ty + transform[2] * tz + transform[3];
          float fy = transform[4] * tx + transform[5] * ty + transform[6] * tz + transform[7];
          float fz = transform[8] * tx + transform[9] * ty + transform[10] * tz + transform[11];
          fwrite(&fx, sizeof(float), 1, fp);
          fwrite(&fy, sizeof(float), 1, fp);
          fwrite(&fz, sizeof(float), 1, fp);

          // Give the points color
          unsigned char r = (unsigned char) (color[0]*255.0f);
          unsigned char g = (unsigned char) (color[1]*255.0f);
          unsigned char b = (unsigned char) (color[2]*255.0f);
          fwrite(&r, sizeof(unsigned char), 1, fp);
          fwrite(&g, sizeof(unsigned char), 1, fp);
          fwrite(&b, sizeof(unsigned char), 1, fp);
        }
      }
    }
  }
  fclose(fp);
}

// Load keypoints, apply transform, and save to point cloud
void keypoints2ply(std::vector<std::vector<float>> &keypoints, const std::string &filename, float* transform, std::vector<float> &color) {

  // Count total number of points in point cloud
  int num_points = keypoints.size()*20;

  // Create header for ply file
  FILE *fp = fopen(filename.c_str(), "w");
  fprintf(fp, "ply\n");
  fprintf(fp, "format binary_little_endian 1.0\n");
  fprintf(fp, "element vertex %d\n", num_points);
  fprintf(fp, "property float x\n");
  fprintf(fp, "property float y\n");
  fprintf(fp, "property float z\n");
  fprintf(fp, "property uchar red\n");
  fprintf(fp, "property uchar green\n");
  fprintf(fp, "property uchar blue\n");
  fprintf(fp, "end_header\n");

  for (int i = 0; i < keypoints.size(); i++) {

    float sx = (float)keypoints[i][0];
    float sy = (float)keypoints[i][1];
    float sz = (float)keypoints[i][2];

    // Apply Rt transform to points
    float fx = transform[0] * sx + transform[1] * sy + transform[2] * sz + transform[3];
    float fy = transform[4] * sx + transform[5] * sy + transform[6] * sz + transform[7];
    float fz = transform[8] * sx + transform[9] * sy + transform[10] * sz + transform[11];
    fwrite(&fx, sizeof(float), 1, fp);
    fwrite(&fy, sizeof(float), 1, fp);
    fwrite(&fz, sizeof(float), 1, fp);

    // Give the points color
    unsigned char r = (unsigned char) (color[0]*255.0f);
    unsigned char g = (unsigned char) (color[1]*255.0f);
    unsigned char b = (unsigned char) (color[2]*255.0f);
    fwrite(&r, sizeof(unsigned char), 1, fp);
    fwrite(&g, sizeof(unsigned char), 1, fp);
    fwrite(&b, sizeof(unsigned char), 1, fp);

    // Draw 3x3 box around keypoint
    for (int kk = -1; kk <= 1; kk++) {
      for (int jj = -1; jj <= 1; jj++) {
        for (int ii = -1; ii <= 1; ii++) {
          int num_edges = 0;
          if (kk == -1 || kk == 1)
            num_edges++;
          if (jj == -1 || jj == 1)
            num_edges++;
          if (ii == -1 || ii == 1)
            num_edges++;
          if (num_edges >= 2) {
            float float_x = (float) (sx + (float)ii*0.01f);
            float float_y = (float) (sy + (float)jj*0.01f);
            float float_z = (float) (sz + (float)kk*0.01f);

            fx = transform[0] * float_x + transform[1] * float_y + transform[2] * float_z + transform[3];
            fy = transform[4] * float_x + transform[5] * float_y + transform[6] * float_z + transform[7];
            fz = transform[8] * float_x + transform[9] * float_y + transform[10] * float_z + transform[11];

            fwrite(&fx, sizeof(float), 1, fp);
            fwrite(&fy, sizeof(float), 1, fp);
            fwrite(&fz, sizeof(float), 1, fp);
            fwrite(&r, sizeof(unsigned char), 1, fp);
            fwrite(&g, sizeof(unsigned char), 1, fp);
            fwrite(&b, sizeof(unsigned char), 1, fp);
          }
        }
      }
    }
  }
  fclose(fp);
}

void get_frag_keypoints(const std::string &frag_dir, std::vector<std::vector<int>> &grid_keypoints,   std::vector<std::vector<float>> &world_keypoints) {

  // If has keypoint file, load it. If not, detect keypoints.
  std::string grid_keypoints_filename = frag_dir + "/keypoints/grid_keypoints.txt";
  std::string world_keypoints_filename = frag_dir + "/keypoints/world_keypoints.txt";
  if (file_exists(grid_keypoints_filename) && file_exists(world_keypoints_filename)) {
    std::cerr << "Loading pre-computed keypoints..." << std::endl;

    // Load keypoints in grid coordinates
    FILE *fp = fopen(grid_keypoints_filename.c_str(), "r");
    int num_keypoints = 0;
    int iret = fscanf(fp, "# number of keypoints: %d", &num_keypoints);
    for (int i = 0; i < num_keypoints; i++) {
      std::vector<int> tmp_grid_keypoint;
      for (int j = 0; j < 3; j++) {
        int tmp_coord;
        iret = fscanf(fp, "%d", &tmp_coord);
        tmp_grid_keypoint.push_back(tmp_coord);
      }
      grid_keypoints.push_back(tmp_grid_keypoint);
    }
    fclose(fp);

    // Load keypoints in world coordinates
    fp = fopen(world_keypoints_filename.c_str(), "r");
    num_keypoints = 0;
    iret = fscanf(fp, "# number of keypoints: %d", &num_keypoints);
    for (int i = 0; i < num_keypoints; i++) {
      std::vector<float> tmp_world_keypoint;
      for (int j = 0; j < 3; j++) {
        float tmp_coord;
        iret = fscanf(fp, "%f", &tmp_coord);
        tmp_world_keypoint.push_back(tmp_coord);
      }
      world_keypoints.push_back(tmp_world_keypoint);
    }
    fclose(fp);

  } else {
    std::cerr << "Keypoints not pre-computed. Detecting keypoints..." << std::endl;
    if (!file_exists(frag_dir + "/keypoints"))
      sys_command("mkdir " + frag_dir + "/keypoints");

    // Load TSDF volume
    std::vector<int> tsdf_dim = load_frag_dim(frag_dir);
    float *tsdf = new float[tsdf_dim[0] * tsdf_dim[1] * tsdf_dim[2]];
    load_frag_tsdf(frag_dir, tsdf, tsdf_dim[0], tsdf_dim[1], tsdf_dim[2]);

    // Detect keypoints in TSDF volume
    float surface_tsdf_threshold = 0.2;
    std::vector<std::vector<int>> detected_keypoints =  detect_keypoints(tsdf, tsdf_dim[0], tsdf_dim[1], tsdf_dim[2], surface_tsdf_threshold, 1.0f, 5, 255.0f);
    grid_keypoints = prune_keypoints(detected_keypoints, 15, tsdf_dim[0], tsdf_dim[1], tsdf_dim[2]);

    // Save keypoints in grid coordinates
    FILE *fp = fopen(grid_keypoints_filename.c_str(), "w");
    fprintf(fp, "# number of keypoints: %d\n", (int) grid_keypoints.size());
    for (int i = 0; i < grid_keypoints.size(); i++)
      fprintf(fp, "%d %d %d\n", grid_keypoints[i][0], grid_keypoints[i][1], grid_keypoints[i][2]);
    fclose(fp);

    // Load grid to world matrix
    float *grid2cam = new float[16];
    float *cam2world = new float[16];
    float *grid2world = new float[16];
    load_frag_ext(frag_dir, grid2cam, cam2world);
    multiply_matrix(cam2world, grid2cam, grid2world);

    // Convert keypoints to world coordinates
    for (int i = 0; i < grid_keypoints.size(); i++) {
      std::vector<float> tmp_world_keypoint;
      float sx = (float) grid_keypoints[i][0];
      float sy = (float) grid_keypoints[i][1];
      float sz = (float) grid_keypoints[i][2];
      float fx = grid2world[0] * sx + grid2world[1] * sy + grid2world[2] * sz + grid2world[3];
      float fy = grid2world[4] * sx + grid2world[5] * sy + grid2world[6] * sz + grid2world[7];
      float fz = grid2world[8] * sx + grid2world[9] * sy + grid2world[10] * sz + grid2world[11];
      tmp_world_keypoint.push_back(fx);
      tmp_world_keypoint.push_back(fy);
      tmp_world_keypoint.push_back(fz);
      world_keypoints.push_back(tmp_world_keypoint);
    }

    // Save keypoints in world coordinates
    fp = fopen(world_keypoints_filename.c_str(), "w");
    fprintf(fp, "# number of keypoints: %d\n", (int) world_keypoints.size());
    for (int i = 0; i < world_keypoints.size(); i++) {
      fprintf(fp, "%.17g %.17g %.17g\n", world_keypoints[i][0], world_keypoints[i][1], world_keypoints[i][2]);
    }
    fclose(fp);

    delete [] tsdf;
    delete [] grid2cam;
    delete [] cam2world;
    delete [] grid2world;

  }
}

void get_frag_features(const std::string &frag_dir, std::vector<std::vector<int>> &grid_keypoints, int feature_size, std::vector<std::vector<float>> &keypoint_features) {


  // If has keypoint features file, load it. If not, find features.
  std::string keypoint_features_filename = frag_dir + "/keypoints/keypoint_features.txt";
  if (file_exists(keypoint_features_filename)) {
    std::cerr << "Loading pre-computed keypoint features..." << std::endl;

    // Load keypoints features
    FILE *fp = fopen(keypoint_features_filename.c_str(), "r");
    int num_keypoints = 0;
    int iret = fscanf(fp, "# number of keypoints: %d", &num_keypoints);
    for (int i = 0; i < num_keypoints; i++) {
      std::vector<float> tmp_keypoint_feature;
      for (int j = 0; j < feature_size; j++) {
        float tmp_val;
        iret = fscanf(fp, "%f", &tmp_val);
        tmp_keypoint_feature.push_back(tmp_val);
      }
      keypoint_features.push_back(tmp_keypoint_feature);
    }
    fclose(fp);

  } else {
    std::cerr << "Keypoint features not pre-computed. Computing features..." << std::endl;
    if (!file_exists(frag_dir + "/keypoints"))
      sys_command("mkdir " + frag_dir + "/keypoints");

    // Load TSDF volume
    std::vector<int> tsdf_dim = load_frag_dim(frag_dir);
    float *tsdf = new float[tsdf_dim[0] * tsdf_dim[1] * tsdf_dim[2]];
    load_frag_tsdf(frag_dir, tsdf, tsdf_dim[0], tsdf_dim[1], tsdf_dim[2]);

    // Find DDD keypoint features
    keypoint_features = ddd_get_keypoint_feat(tsdf, tsdf_dim[0], tsdf_dim[1], tsdf_dim[2], grid_keypoints, 15, true);

    // Save keypoints features
    FILE *fp = fopen(keypoint_features_filename.c_str(), "w");
    fprintf(fp, "# number of keypoints: %d\n", (int) keypoint_features.size());
    for (int i = 0; i < keypoint_features.size(); i++) {
      for (int j = 0; j < keypoint_features[i].size(); j++)
        fprintf(fp, "%.17g ", keypoint_features[i][j]);
      fprintf(fp, "\n");
    }
    fclose(fp);

    delete [] tsdf;

  }
}

void compare_frag_features(std::vector<std::vector<float>> &keypoint_features1, std::vector<std::vector<float>> &keypoint_features2, const std::string &frag_cache_dir, const std::string &frag_pair_name, std::vector<std::vector<float>> &score_matrix1, std::vector<std::vector<float>> &score_matrix2) {

  // If cache directory doesn't exist, create it
  if (!file_exists(frag_cache_dir))
    sys_command("mkdir " + frag_cache_dir);

  std::string score_matrix1_filename = frag_cache_dir + "/" + frag_pair_name + "_scores1.txt";
  std::string score_matrix2_filename = frag_cache_dir + "/" + frag_pair_name + "_scores2.txt";
  // Check to see if cache contains score matrices for the comparison
  if (file_exists(score_matrix1_filename) && file_exists(score_matrix2_filename)) {

    // Load score matrix 1
    FILE *fp = fopen(score_matrix1_filename.c_str(), "r");
    int num_rows = 0;
    int num_cols = 0;
    int iret = fscanf(fp, "# number of keypoint comparisons: %d x %d", &num_rows, &num_cols);
    for (int i = 0; i < num_rows; i++) {
      std::vector<float> tmp_scores;
      for (int j = 0; j < num_cols; j++) {
        float tmp_val;
        iret = fscanf(fp, "%f", &tmp_val);
        tmp_scores.push_back(tmp_val);
      }
      score_matrix1.push_back(tmp_scores);
    }
    fclose(fp);

    // Load score matrix 2
    fp = fopen(score_matrix2_filename.c_str(), "r");
    num_rows = 0;
    num_cols = 0;
    iret = fscanf(fp, "# number of keypoint comparisons: %d x %d", &num_rows, &num_cols);
    for (int i = 0; i < num_rows; i++) {
      std::vector<float> tmp_scores;
      for (int j = 0; j < num_cols; j++) {
        float tmp_val;
        iret = fscanf(fp, "%f", &tmp_val);
        tmp_scores.push_back(tmp_val);
      }
      score_matrix2.push_back(tmp_scores);
    }
    fclose(fp);

  } else {

    ddd_compare_feat(keypoint_features1, keypoint_features2, score_matrix1, true);
    ddd_compare_feat(keypoint_features2, keypoint_features1, score_matrix2, true);

    // Save score matrix 1
    FILE *fp = fopen(score_matrix1_filename.c_str(), "w");
    fprintf(fp, "# number of keypoint comparisons: %d x %d\n", (int) score_matrix1.size(), (int) score_matrix2.size());
    for (int i = 0; i < score_matrix1.size(); i++) {
      for (int j = 0; j < score_matrix1[i].size(); j++)
        fprintf(fp, "%.17g ", score_matrix1[i][j]);
      fprintf(fp, "\n");
    }
    fclose(fp);

    // Save score matrix 2
    fp = fopen(score_matrix2_filename.c_str(), "w");
    fprintf(fp, "# number of keypoint comparisons: %d x %d\n", (int) score_matrix2.size(), (int) score_matrix1.size());
    for (int i = 0; i < score_matrix2.size(); i++) {
      for (int j = 0; j < score_matrix2[i].size(); j++)
        fprintf(fp, "%.17g ", score_matrix2[i][j]);
      fprintf(fp, "\n");
    }
    fclose(fp);

  }
}

void save_frag_world_keypoints_cache(std::vector<std::vector<float>> &world_keypoints1, std::vector<std::vector<float>> &world_keypoints2, const std::string &frag_cache_dir, const std::string &frag_pair_name) {

  std::string cache_world_keypoints_filename_1 = frag_cache_dir + "/" + frag_pair_name + "_keypoints1.txt";
  std::string cache_world_keypoints_filename_2 = frag_cache_dir + "/" + frag_pair_name + "_keypoints2.txt";

  if (!file_exists(cache_world_keypoints_filename_1) || !file_exists(cache_world_keypoints_filename_2)) {

    // Save keypoints from first fragment to cache
    FILE *fp = fopen(cache_world_keypoints_filename_1.c_str(), "w");
    fprintf(fp, "# number of keypoints: %d\n", (int) world_keypoints1.size());
    for (int i = 0; i < world_keypoints1.size(); i++) {
      fprintf(fp, "%.17g %.17g %.17g\n", world_keypoints1[i][0], world_keypoints1[i][1], world_keypoints1[i][2]);
    }
    fclose(fp);

    // Save keypoints from second fragment to cache
    fp = fopen(cache_world_keypoints_filename_2.c_str(), "w");
    fprintf(fp, "# number of keypoints: %d\n", (int) world_keypoints2.size());
    for (int i = 0; i < world_keypoints2.size(); i++) {
      fprintf(fp, "%.17g %.17g %.17g\n", world_keypoints2[i][0], world_keypoints2[i][1], world_keypoints2[i][2]);
    }
    fclose(fp);

  }
}