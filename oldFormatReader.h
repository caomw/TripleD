

void checkout_tsdf(const std::string &scene_name, float* tsdf_volume, int x1, int x2, int y1, int y2, int z1, int z2) {
  // Reads in volume data (xyz order float format) specified by indicies (inclusive)
  // TSDF volume file format:
  //     start index of volume 1*(int32)
  //     number of elements (N) in volume 1*(int32)
  //     volume values N*(float)

  std::string filename = scene_name + ".tsdf";

  std::ifstream inFile(filename, std::ios::binary | std::ios::in);
  int offset;
  int num_elements;
  inFile.read((char*)&offset, sizeof(int));
  inFile.read((char*)&num_elements, sizeof(int));
  int end_idx = offset + num_elements - 1;

  int x_dim = x2 - x1 + 1;
  int y_dim = y2 - y1 + 1;
  // int z_dim = z2 - z1 + 1;

  // for (int z = z1; z <= z2; z++) {
  //   for (int y = y1; y <= y2; y++) {
  //     for (int x = x1; x <= x2; x++) {
  //       int volume_idx = z * 512 * 512 + y * 512 + x;
  //       if (volume_idx > end_idx || volume_idx < offset)
  //         tsdf_volume[(z - z1)*y_dim * x_dim + (y - y1)*x_dim + (x - x1)] = 1.0f;
  //       else {
  //         inFile.seekg(8 + sizeof(float) * (volume_idx - offset));
  //         // float tsdf_value;
  //         inFile.read((char*)&tsdf_volume[(z - z1)*y_dim * x_dim + (y - y1)*x_dim + (x - x1)], sizeof(float));
  //         // tsdf_volume[(z-z1)*y_dim*x_dim+(y-y1)*x_dim+(x-x1)] = tsdf_value;
  //       }
  //     }
  //   }
  // }

  for (int z = z1; z <= z2; z++) {
    for (int y = y1; y <= y2; y++) {
      int volume_idx = z * 512 * 512 + y * 512 + x1;
      // If can't load entire x row
      if ((volume_idx - x1 + x2) > end_idx || volume_idx < offset) {
        for (int x = x1; x <= x2; x++) {
          volume_idx = z * 512 * 512 + y * 512 + x;
          if (volume_idx > end_idx || volume_idx < offset)
            tsdf_volume[(z - z1)*y_dim * x_dim + (y - y1)*x_dim + (x - x1)] = 1.0f;
          else {
            inFile.seekg(8 + sizeof(float) * (volume_idx - offset));
            // float tsdf_value;
            inFile.read((char*)&tsdf_volume[(z - z1)*y_dim * x_dim + (y - y1)*x_dim + (x - x1)], sizeof(float));
            // tsdf_volume[(z-z1)*y_dim*x_dim+(y-y1)*x_dim+(x-x1)] = tsdf_value;
          }
        }
        // If can load entire x row
      } else {
        inFile.seekg(8 + sizeof(float) * (volume_idx - offset));
        inFile.read((char*)&tsdf_volume[(z - z1)*y_dim * x_dim + (y - y1)*x_dim], x_dim * sizeof(float));
      }
    }
  }

  inFile.close();

}

void checkout_ext(const std::string &scene_name, float* ext_mat) {

  std::string filename = scene_name + "_ext.txt";

  int iret;
  FILE *fp = fopen(filename.c_str(), "r");
  for (int i = 0; i < 16; i++)
    iret = fscanf(fp, "%f", &ext_mat[i]);
  fclose(fp);
}

////////////////////////////////////////////////////////////////////////////////

void checkout_keypts(const std::string &scene_name, std::vector<std::vector<float>> &keypoints) {

  std::string filename = scene_name + "_pts.txt";

  FILE *fp = fopen(filename.c_str(), "r");
  int iret;
  float num_points;
  iret = fscanf(fp, "%f", &num_points);
  for (int i = 0; i < num_points; i++) {
    std::vector<float> tmp_keypt;
    tmp_keypt.resize(3);
    iret = fscanf(fp, "%f", &tmp_keypt[0]);
    iret = fscanf(fp, "%f", &tmp_keypt[1]);
    iret = fscanf(fp, "%f", &tmp_keypt[2]);
    keypoints.push_back(tmp_keypt);
  }
  fclose(fp);
}

void convert_tsdf_old_to_new(const std::string &src, const std::string &dst) {

  std::vector<std::vector<float>> grid_keypoints;
  checkout_keypts(src, grid_keypoints);
  if (grid_keypoints.size() == 0)
    return;

  // Read TSDF from old format
  float *tsdf = new float[512 * 512 * 1024];
  checkout_tsdf(src, tsdf, 0, 512 - 1, 0, 512 - 1, 0, 1024 - 1);

  // Find dimensions of TSDF volume
  int x_min = 511, x_max = 0, y_min = 511, y_max = 0, z_min = 1023, z_max = 0; // inclusive bounds
  for (int z = 0; z < 1024; z++)
    for (int y = 0; y < 512; y++)
      for (int x = 0; x < 512; x++) {
        int i = z * 512 * 512 + y * 512 + x;
        if (tsdf[i] < 1) {
          x_max = std::max(x, x_max);
          y_max = std::max(y, y_max);
          z_max = std::max(z, z_max);
        }
      }
  for (int z = 1023; z >= 0; z--)
    for (int y = 511; y >= 0; y--)
      for (int x = 511; x >= 0; x--) {
        int i = z * 512 * 512 + y * 512 + x;
        if (tsdf[i] < 1) {
          x_min = std::min(x, x_min);
          y_min = std::min(y, y_min);
          z_min = std::min(z, z_min);
        }
      }
  int x_dim = x_max - x_min + 1;
  int y_dim = y_max - y_min + 1;
  int z_dim = z_max - z_min + 1;

  // // Debug
  // std::cout << x_min << std::endl;
  // std::cout << x_max << std::endl;
  // std::cout << y_min << std::endl;
  // std::cout << y_max << std::endl;
  // std::cout << z_min << std::endl;
  // std::cout << z_max << std::endl;

  // Create new smaller TSDF volume
  float *new_tsdf = new float[x_dim * y_dim * z_dim];
  for (int z = 0; z < z_dim; z++)
    for (int y = 0; y < y_dim; y++)
      for (int x = 0; x < x_dim; x++)
        new_tsdf[z * y_dim * x_dim + y * x_dim + x] = tsdf[(z + z_min) * 512 * 512 + (y + y_min) * 512 + (x + x_min)];

  sys_command("mkdir " + dst);

  save_frag_tsdf(dst, new_tsdf, x_dim, y_dim, z_dim); // function found in fragments/io.h

  // float *cam2world = new float[16];
  // checkout_ext(src, cam2world);

  float *grid2cam = new float[16];
  grid2cam[0] = 0.01f;  grid2cam[1] = 0.0f;   grid2cam[2] = 0.0f;    grid2cam[3] = -2.55f + ((float)x_min)*0.01f; 
  grid2cam[4] = 0.0f;   grid2cam[5] = 0.01f;  grid2cam[6] = 0.0f;    grid2cam[7] = -2.55f + ((float)y_min)*0.01f; 
  grid2cam[8] = 0.0f;   grid2cam[9] = 0.0f;   grid2cam[10] = 0.01f;  grid2cam[11] = -0.49f + ((float)z_min)*0.01f; 
  grid2cam[12] = 0.0f;  grid2cam[13] = 0.0f;  grid2cam[14] = 0.0f;   grid2cam[15] = 1.0f; 

  // float *grid2firstcam = new float[16];
  // multiply_matrix(cam2world, grid2cam, grid2firstcam);

  // grid2world[0] = 1.0f;
  // grid2world[1] = 0.0f;
  // grid2world[2] = 0.0f;
  // grid2world[3] = 0.0f;
  // grid2world[4] = 0.0f;
  // grid2world[5] = 1.0f;
  // grid2world[6] = 0.0f;
  // grid2world[7] = 0.0f;
  // grid2world[8] = 0.0f;
  // grid2world[9] = 0.0f;
  // grid2world[10] = 1.0f;
  // grid2world[11] = 0.0f;
  // grid2world[12] = 0.0f;
  // grid2world[13] = 0.0f;
  // grid2world[14] = 0.0f;
  // grid2world[15] = 1.0f;

  // std::string grid2cam_filename = dst + "/grid2cam.txt";
  // FILE *fp = fopen(grid2cam_filename.c_str(), "w");
  // fprintf(fp, "%.17g %.17g %.17g %.17g\n", grid2cam[0], grid2cam[1], grid2cam[2], grid2cam[3]);
  // fprintf(fp, "%.17g %.17g %.17g %.17g\n", grid2cam[4], grid2cam[5], grid2cam[6], grid2cam[7]);
  // fprintf(fp, "%.17g %.17g %.17g %.17g\n", grid2cam[8], grid2cam[9], grid2cam[10], grid2cam[11]);
  // fprintf(fp, "%.17g %.17g %.17g %.17g\n", grid2cam[12], grid2cam[13], grid2cam[14], grid2cam[15]);
  // fclose(fp);

  // std::string cam2world_filename = dst + "/cam2world.txt";
  // fp = fopen(cam2world_filename.c_str(), "w");
  // fprintf(fp, "%.17g %.17g %.17g %.17g\n", cam2world[0], cam2world[1], cam2world[2], cam2world[3]);
  // fprintf(fp, "%.17g %.17g %.17g %.17g\n", cam2world[4], cam2world[5], cam2world[6], cam2world[7]);
  // fprintf(fp, "%.17g %.17g %.17g %.17g\n", cam2world[8], cam2world[9], cam2world[10], cam2world[11]);
  // fprintf(fp, "%.17g %.17g %.17g %.17g\n", cam2world[12], cam2world[13], cam2world[14], cam2world[15]);
  // fclose(fp);

  std::string grid2cam_filename = dst + "/grid2cam.txt";
  FILE *fp = fopen(grid2cam_filename.c_str(), "w");
  fprintf(fp, "%.17g %.17g %.17g %.17g\n", grid2cam[0], grid2cam[1], grid2cam[2], grid2cam[3]);
  fprintf(fp, "%.17g %.17g %.17g %.17g\n", grid2cam[4], grid2cam[5], grid2cam[6], grid2cam[7]);
  fprintf(fp, "%.17g %.17g %.17g %.17g\n", grid2cam[8], grid2cam[9], grid2cam[10], grid2cam[11]);
  fprintf(fp, "%.17g %.17g %.17g %.17g\n", grid2cam[12], grid2cam[13], grid2cam[14], grid2cam[15]);
  fclose(fp);

  tsdf2ply(dst + "/pointcloud.ply", new_tsdf, 0.2, x_dim, y_dim, z_dim, grid2cam);

  std::vector<std::vector<int>> grid_keypoints1;
  std::vector<std::vector<float>> cam_keypoints1;
  get_frag_keypoints(dst, grid_keypoints1, cam_keypoints1);

  float* transform = new float[16];
  load_identity_matrix(transform);
  std::vector<float> color_red;
  color_red.push_back(255.0f / 255.0f);
  color_red.push_back(0.0f / 255.0f);
  color_red.push_back(0.0f / 255.0f);
  keypoints2ply(cam_keypoints1, dst + "/keypoints/harris/keypoints.ply", transform, color_red);
  
  std::vector<std::vector<float>> keypoint_features1;
  get_frag_features(dst, grid_keypoints1, 2048, keypoint_features1);

  // std::vector<std::vector<float>> grid_keypoints;
  // checkout_keypts(src, grid_keypoints);
  // for (int i = 0; i < grid_keypoints.size(); i++) {
  //   grid_keypoints[i][0]  = grid_keypoints[i][0] - x_min;
  //   grid_keypoints[i][1]  = grid_keypoints[i][1] - y_min;
  //   grid_keypoints[i][2]  = grid_keypoints[i][2] - z_min;
  //   std::cout << grid_keypoints[i][0] << " " << grid_keypoints[i][1] << " " << grid_keypoints[i][2] << std::endl;
  // }


  // std::vector<std::vector<float>> cam_keypoints;


  // if (!file_exists(dst + "/keypoints"))
  //     sys_command("mkdir " + dst + "/keypoints");

  // std::string grid_keypoints_filename = dst + "/keypoints/grid_keypoints.txt";
  // std::string cam_keypoints_filename = dst + "/keypoints/cam_keypoints.txt";

  // // Save keypoints in grid coordinates
  // fp = fopen(grid_keypoints_filename.c_str(), "w");
  // fprintf(fp, "# number of keypoints: %d\n", (int) grid_keypoints.size());
  // for (int i = 0; i < grid_keypoints.size(); i++)
  //   fprintf(fp, "%d %d %d\n", (int)grid_keypoints[i][0], (int)grid_keypoints[i][1], (int)grid_keypoints[i][2]);
  // fclose(fp);

  // // Convert keypoints to world coordinates
  // for (int i = 0; i < grid_keypoints.size(); i++) {
  //   std::vector<float> tmp_cam_keypoint;
  //   float sx = (float) grid_keypoints[i][0];
  //   float sy = (float) grid_keypoints[i][1];
  //   float sz = (float) grid_keypoints[i][2];
  //   float fx = grid2world[0] * sx + grid2world[1] * sy + grid2world[2] * sz + grid2world[3];
  //   float fy = grid2world[4] * sx + grid2world[5] * sy + grid2world[6] * sz + grid2world[7];
  //   float fz = grid2world[8] * sx + grid2world[9] * sy + grid2world[10] * sz + grid2world[11];
  //   tmp_cam_keypoint.push_back(fx);
  //   tmp_cam_keypoint.push_back(fy);
  //   tmp_cam_keypoint.push_back(fz);
  //   cam_keypoints.push_back(tmp_cam_keypoint);
  // }

  // // Save keypoints in world coordinates
  // fp = fopen(cam_keypoints_filename.c_str(), "w");
  // fprintf(fp, "# number of keypoints: %d\n", (int) cam_keypoints.size());
  // for (int i = 0; i < cam_keypoints.size(); i++) {
  //   fprintf(fp, "%.17g %.17g %.17g\n", cam_keypoints[i][0], cam_keypoints[i][1], cam_keypoints[i][2]);
  // }
  // fclose(fp);



  // }


  delete [] tsdf;
  delete [] new_tsdf;

}