#include <uWS/uWS.h>
#include <iostream>
#include "json.hpp"
#include <math.h>
#include "particle_filter.h"
#include "helper_functions.h"

using namespace std;

#include <random>

int testInit() {
  ParticleFilter pf;
  double zero_sigma_pos [3] = {0., 0., 0.};
  const double gps_x = 4983;
  const double gps_y = 5029;
  const double theta = 1.201;
  const int numberOfParticles=100;
  pf.setNumberOfParticles(numberOfParticles);
  
  pf.init(gps_x, gps_y, theta, zero_sigma_pos);
  assert(pf.particles.size()==numberOfParticles);
  for (int p=0; p<pf.particles.size(); p++) {
    assert(areSame(pf.particles[p].x, gps_x));
    assert(areSame(pf.particles[p].y, gps_y));
    assert(areSame(pf.particles[p].theta, theta));
  }
  pf.deleteAllParticlesInFilter();
  double sigma_pos [3] = {0.1, 0.1, 0.1};
  pf.init(gps_x, gps_y, theta, sigma_pos);
  assert(pf.particles.size()==numberOfParticles);
  for (int p=0; p<pf.particles.size(); p++) {
    assert(!areSame(pf.particles[p].x, gps_x));
    assert(!areSame(pf.particles[p].y, gps_y));
    assert(!areSame(pf.particles[p].theta, theta));
    const double epsilon=0.5;
    assert(areSame(pf.particles[p].x, gps_x, epsilon));
    assert(areSame(pf.particles[p].y, gps_y, epsilon));
    assert(areSame(pf.particles[p].theta, theta, epsilon));
  }
  
  assert(pf.initialized());
  
  return 0;
}

int testPredict() {
  
  ParticleFilter pf0;
  
  const int numberOfParticles=10;
  pf0.setNumberOfParticles(numberOfParticles);
  
  const double init_x = 102;
  const double init_y = 65;
  const double init_theta = (5.*M_PI)/8.;
  const double zero_sigma_pos [3] = {0., 0., 0.};
  
  pf0.init(init_x, init_y, init_theta, zero_sigma_pos);
  assert(pf0.particles.size()==numberOfParticles);
  for (int p=0; p<pf0.particles.size(); p++) {
    assert(areSame(pf0.particles[p].x, init_x));
    assert(areSame(pf0.particles[p].y, init_y));
    assert(areSame(pf0.particles[p].theta, init_theta));
  }
  
  pf0.prediction(0./*deltaT*/, zero_sigma_pos ,  1., 1.);
  for (int p=0; p<pf0.particles.size(); p++) {
    assert(areSame(pf0.particles[p].x, init_x));
    assert(areSame(pf0.particles[p].y, init_y));
    assert(areSame(pf0.particles[p].theta, init_theta));
  }
  const double calculated_x = 97.592046082729055;
  const double calculated_y = 75.07741997215382;
  const double calculated_theta = (51.*M_PI)/80.;
  
  pf0.prediction(0.1/*deltaT*/, zero_sigma_pos ,  110., M_PI/8.);
  for (int p=0; p<pf0.particles.size(); p++) {
    assert(areSame(pf0.particles[p].x, calculated_x));
    assert(areSame(pf0.particles[p].y, calculated_y));
    assert(areSame(pf0.particles[p].theta, calculated_theta));
  }
  pf0.prediction(0./*deltaT*/, (double[]) {0.1, 0.1, 0.1} ,  1., 1.);
  // deltaT==0. -> x,y,theta can't change
  for (int p=0; p<pf0.particles.size(); p++) {
    assert(!areSame(pf0.particles[p].x, calculated_x));
    assert(!areSame(pf0.particles[p].y, calculated_y));
    assert(!areSame(pf0.particles[p].theta, calculated_theta));
    const double epsilon=0.5;
    assert(areSame(pf0.particles[p].x, calculated_x, epsilon));
    assert(areSame(pf0.particles[p].y, calculated_y, epsilon));
    assert(areSame(pf0.particles[p].theta, calculated_theta, epsilon));
  }
  
  // testing 0 yawRate
  
  const double deltaT = 1.0;
  const double velocity = 100.0;
  const double initThetaXOnly = 0.;// particles all facing right
  const double zeroYawRate = 0.;// particles only move in facing direction
  
  ParticleFilter pf1;
  pf1.setNumberOfParticles(numberOfParticles);
  
  pf1.init(init_x, init_y, initThetaXOnly, zero_sigma_pos);
  assert(pf1.particles.size()==numberOfParticles);
  for (int p=0; p<pf1.particles.size(); p++) {
    assert(areSame(pf1.particles[p].x, init_x));
    assert(areSame(pf1.particles[p].y, init_y));
    assert(areSame(pf1.particles[p].theta, initThetaXOnly));
  }
  pf1.prediction(deltaT, zero_sigma_pos ,  velocity, zeroYawRate);
  for (int p=0; p<pf1.particles.size(); p++) {
    assert(areSame(pf1.particles[p].x, init_x+(velocity*deltaT)));// v*dt in x only
    assert(areSame(pf1.particles[p].y, init_y));
    assert(areSame(pf1.particles[p].theta, initThetaXOnly));
  }
  
  ParticleFilter pf2;
  pf2.setNumberOfParticles(numberOfParticles);
  const double initThetaYOnly = M_PI/2.;// particles all facing up
  
  pf2.init(init_x, init_y, initThetaYOnly, zero_sigma_pos);
  assert(pf2.particles.size()==numberOfParticles);
  for (int p=0; p<pf2.particles.size(); p++) {
    assert(areSame(pf2.particles[p].x, init_x));
    assert(areSame(pf2.particles[p].y, init_y));
    assert(areSame(pf2.particles[p].theta, initThetaYOnly));
  }
  pf2.prediction(deltaT, zero_sigma_pos ,  velocity, zeroYawRate);
  for (int p=0; p<pf2.particles.size(); p++) {
    assert(areSame(pf2.particles[p].x, init_x));
    assert(areSame(pf2.particles[p].y, init_y)+(velocity*deltaT));// v*dt in y only
    assert(areSame(pf2.particles[p].theta, initThetaYOnly));
  }
  
  ParticleFilter pf3;
  pf3.setNumberOfParticles(2);
  
  pf3.init(init_x, init_y, init_theta, zero_sigma_pos);
  assert(pf3.particles.size()==2);
  for (int p=0; p<pf3.particles.size(); p++) {
    assert(areSame(pf3.particles[p].x, init_x));
    assert(areSame(pf3.particles[p].y, init_y));
    assert(areSame(pf3.particles[p].theta, init_theta));
  }
  
  // noise only
  double sigma_pos [3] = {0.3, 0.3, 0.01};
  
  pf3.prediction(0.0/* deltaT */, sigma_pos ,  0./* velocity */, 0./* yaw rate */);
  for (int p=0; p<pf3.particles.size(); p++) {
    assert(areSame(pf3.particles[p].x, init_x, 3.*sigma_pos[0]));
    assert(areSame(pf3.particles[p].y, init_y, 3.*sigma_pos[1]));
    assert(areSame(pf3.particles[p].theta, init_theta, 3.*sigma_pos[2]));
  }
  
  return 0;
}

int testError() {
  // (5.2 m, 19.3 m, pi/16
  const double gt_x=5.2;
  const double gt_y=19.3;
  const double gt_theta=M_PI/16.;
  // 5 m, 18.7 m, pi/8
  const double pf_x=5.;
  const double pf_y=18.7;
  const double pf_theta=M_PI/8.;
  
  const double *error = getError(gt_x, gt_y, gt_theta, pf_x, pf_y, pf_theta);
  
  const double positionError = rsmePosition(error);
  assert(areSame(positionError, 0.63245553203367733));
  
  const double angleError = rsmeTheta(error);
  assert(areSame(angleError, M_PI/16.));
  
  return 0;
}

int testTransform() {
  
  Particle particle;
  particle.x=4.;
  particle.y=5.;
  particle.theta = -M_PI/2.;
  
  LandmarkObs observation1;
  observation1.x=2.;
  observation1.y=2.;
  LandmarkObs transformedObservation1=ParticleFilter::transform(particle, observation1);
  assert (areSame(transformedObservation1.x, 6.));
  assert (areSame(transformedObservation1.y, 3.));
  
  LandmarkObs observation2;
  observation2.x=3.;
  observation2.y=-2.;
  LandmarkObs transformedObservation2=ParticleFilter::transform(particle, observation2);
  assert (areSame(transformedObservation2.x, 2.));
  assert (areSame(transformedObservation2.y, 2.));
  
  LandmarkObs observation3;
  observation3.x=0.;
  observation3.y=-4.;
  LandmarkObs transformedObservation3=ParticleFilter::transform(particle, observation3);
  assert (areSame(transformedObservation3.x, 0.));
  assert (areSame(transformedObservation3.y, 5.));
  
  std::vector<LandmarkObs> observations;
  
  observations.push_back(observation1);
  observations.push_back(observation2);
  observations.push_back(observation3);
  
  std::vector<LandmarkObs> transformedObservations = ParticleFilter::transform(particle, observations);
  
  assert (areSame(transformedObservations[0].x, 6.));
  assert (areSame(transformedObservations[0].y, 3.));
  
  assert (areSame(transformedObservations[1].x, 2.));
  assert (areSame(transformedObservations[1].y, 2.));
  
  assert (areSame(transformedObservations[2].x, 0.));
  assert (areSame(transformedObservations[2].y, 5.));
  
  return 0;
}

int testDataAssociation() {
  
  std::vector<LandmarkObs> landmarks;
  
  LandmarkObs landmark1;
  landmark1.id=1;
  landmark1.x=5.;
  landmark1.y=3.;
  landmarks.push_back(landmark1);
  
  LandmarkObs landmark2;
  landmark2.id=2;
  landmark2.x=2.;
  landmark2.y=1.;
  landmarks.push_back(landmark2);
  
  LandmarkObs landmark3;
  landmark3.id=3;
  landmark3.x=6.;
  landmark3.y=1.;
  landmarks.push_back(landmark3);
  
  LandmarkObs landmark4;
  landmark4.id=4;
  landmark4.x=7.;
  landmark4.y=4.;
  landmarks.push_back(landmark4);
  
  LandmarkObs landmark5;
  landmark5.id=5;
  landmark5.x=4.;
  landmark5.y=7.;
  landmarks.push_back(landmark5);
  
  std::vector<LandmarkObs> observations;
  
  LandmarkObs observation1;
  observation1.x=2.;
  observation1.y=2.;
  observations.push_back(observation1);
  
  LandmarkObs observation2;
  observation2.x=3.;
  observation2.y=-2.;
  observations.push_back(observation2);
  
  LandmarkObs observation3;
  observation3.x=0.;
  observation3.y=-4.;
  observations.push_back(observation3);
  
  Particle particle;
  particle.x=4.;
  particle.y=5.;
  particle.theta = -M_PI/2.;
  
  ParticleFilter pf;
  std::vector<LandmarkObs> /* xMap, yMap */  transformedObservations=ParticleFilter::transform(particle, observations);
  pf.dataAssociation(transformedObservations, landmarks);
  
  assert(transformedObservations[0].id==0);// 0 -> L1
  assert(transformedObservations[1].id==1);// 1 -> L2
  assert(transformedObservations[2].id==1);// 1 -> L2
  
  return 0;
}

int testParticleWeight() {
  
  std::vector<LandmarkObs> landmarks;
  
  LandmarkObs landmark1;
  landmark1.id=1;
  landmark1.x=5.;
  landmark1.y=3.;
  landmarks.push_back(landmark1);
  
  LandmarkObs landmark2;
  landmark2.id=2;
  landmark2.x=2.;
  landmark2.y=1.;
  landmarks.push_back(landmark2);
  
  LandmarkObs landmark3;
  landmark3.id=3;
  landmark3.x=6.;
  landmark3.y=1.;
  landmarks.push_back(landmark3);
  
  LandmarkObs landmark4;
  landmark4.id=4;
  landmark4.x=7.;
  landmark4.y=4.;
  landmarks.push_back(landmark4);
  
  LandmarkObs landmark5;
  landmark5.id=5;
  landmark5.x=4.;
  landmark5.y=7.;
  landmarks.push_back(landmark5);
  
  std::vector<LandmarkObs> observations;
  
  LandmarkObs observation1;
  observation1.x=2.;
  observation1.y=2.;
  observations.push_back(observation1);
  
  LandmarkObs observation2;
  observation2.x=3.;
  observation2.y=-2.;
  observations.push_back(observation2);
  
  LandmarkObs observation3;
  observation3.x=0.;
  observation3.y=-4.;
  observations.push_back(observation3);
  
  Particle particle;
  particle.x=4.;
  particle.y=5.;
  particle.theta = -M_PI/2.;
  
  std::vector<LandmarkObs> /* xMap, yMap */  transformedObservations=ParticleFilter::transform(particle, observations);
  for (int tO=0; tO<transformedObservations.size(); tO++) {
    assert(transformedObservations[tO].id<0); // assert transformedObservations are unlinked
  }
  const double sigmas[] = {0.3, 0.3};
  
  //   static const double particleWeight(const double theMeanX, const double theMeanY, const Particle theParticle, const double theParticleSigmas[]/* sx, sy */) {
  
  const double particleWeight1 = ParticleFilter::particleWeight(transformedObservations[0] /* observation 1*/, landmark1, sigmas);
  assert(areSame(particleWeight1, 6.84e-3));
  
  const double particleWeight2 = ParticleFilter::particleWeight(transformedObservations[1] /* observation 2*/, landmark2, sigmas);
  assert(areSame(particleWeight2, 6.84e-3));
  
  const double particleWeight3 = ParticleFilter::particleWeight(transformedObservations[2] /* observation 3*/, landmark2, sigmas);
  assert(areSame(particleWeight3, 9.83e-49));
  
  ParticleFilter pf;
  pf.dataAssociation(transformedObservations, landmarks);
  for (int tO=0; tO<transformedObservations.size(); tO++) {
    assert(transformedObservations[tO].id>=0); // assert transformedObservations are linked, links are tested in testDataAssociation
  }
  const double particleWeight = ParticleFilter::particleWeight(transformedObservations, landmarks, sigmas);
  assert(areSame(particleWeight, 4.60E-53));
  
  return 0;
}

int testFilterLandmarkMap() {
  
  std::vector<Map::single_landmark_s> landmarks;// in map coordinates
  
  Map::single_landmark_s single_landmark1;
  single_landmark1.id_i=1;
  single_landmark1.x_f=5.;// 1.
  single_landmark1.y_f=3.;// 2.
  landmarks.push_back(single_landmark1);
  
  Map::single_landmark_s single_landmark2;
  single_landmark2.id_i=2;
  single_landmark2.x_f=2.;// 2.
  single_landmark2.y_f=1.;// 4.
  landmarks.push_back(single_landmark2);
  
  Map::single_landmark_s single_landmark3;
  single_landmark3.id_i=3;
  single_landmark3.x_f=6.;// 2.
  single_landmark3.y_f=1.;// 4.
  landmarks.push_back(single_landmark3);
  
  Map::single_landmark_s single_landmark4;
  single_landmark4.id_i=4;
  single_landmark4.x_f=7.;// 3.
  single_landmark4.y_f=4.;// 1.
  landmarks.push_back(single_landmark4);
  
  Map::single_landmark_s single_landmark5;
  single_landmark5.id_i=5;
  single_landmark5.x_f=4.;// 0.
  single_landmark5.y_f=7.;// 2.
  landmarks.push_back(single_landmark5);
  
  Map map;
  map.landmark_list=landmarks;
  
  Particle particle;// in map coordinates
  particle.x=4.;
  particle.y=5.;
  particle.theta = -M_PI/2.;
  
  std::vector<LandmarkObs> allLandmarks=ParticleFilter::filterLandmarkMap(map, particle, 5.); // 5 should keep all landmarks
  assert(allLandmarks.size()==map.landmark_list.size() && allLandmarks.size()>0);
  
  std::vector<LandmarkObs> someLandmarks=ParticleFilter::filterLandmarkMap(map, particle, 3.); // 3 should keep 1 & 5 landmarks
  assert(someLandmarks.size()==2 && someLandmarks[0].id==1 && someLandmarks[1].id==5);
  
  return 0;
}

int testResample() {
  
  ParticleFilter pf;
  double zero_sigma_pos [3] = {0., 0., 0.};
  const double gps_x = 4983;
  const double gps_y = 5029;
  const double theta = 1.201;
  
  const int numberOfParticles=10;
  pf.setNumberOfParticles(numberOfParticles);
  const int firstParticleIdAtTestStart = pf.currentParticleId();
  
  pf.init(gps_x, gps_y, theta, zero_sigma_pos);
  assert(pf.particles.size()==numberOfParticles);
  for (int p=0; p<pf.particles.size(); p++) {
    assert(areSame(pf.particles[p].x, gps_x));
    assert(areSame(pf.particles[p].y, gps_y));
    assert(areSame(pf.particles[p].theta, theta));
    assert(areSame(pf.particles[p].weight, 1.));
  }
  
  const int lastParticleIdAtTestStart = pf.currentParticleId();
  
  pf.resample();
  // pf.particles probably has Particle repeats now
  assert(pf.particles.size()==numberOfParticles);
  for (int p=0; p<pf.particles.size(); p++) {
    assert(areSame(pf.particles[p].x, gps_x));
    assert(areSame(pf.particles[p].y, gps_y));
    assert(areSame(pf.particles[p].theta, theta));
    assert(areSame(pf.particles[p].weight, 1.));
  }
  
  const int numberOfParticlesToRemove=5;
  assert(numberOfParticlesToRemove<numberOfParticles);
  for (int p=numberOfParticlesToRemove; p<pf.particles.size(); p++) {// flag
    pf.particles[p].weight=0.;
  }
  for (int p=0; p<numberOfParticlesToRemove; p++) {
    assert(areSame(pf.particles[p].x, gps_x));
    assert(areSame(pf.particles[p].y, gps_y));
    assert(areSame(pf.particles[p].theta, theta));
    assert(areSame(pf.particles[p].weight, 1.));
  }
  for (int p=numberOfParticlesToRemove; p<numberOfParticles; p++) {
    assert(areSame(pf.particles[p].x, gps_x));
    assert(areSame(pf.particles[p].y, gps_y));
    assert(areSame(pf.particles[p].theta, theta));
    assert(areSame(pf.particles[p].weight, 0.));
  }
  
  pf.resample();
  
  assert(pf.particles.size()==numberOfParticles);
  for (int p=0; p<pf.particles.size(); p++) {
    const Particle& particle=pf.particles[p];
    assert(areSame(pf.particles[p].x, gps_x));
    assert(areSame(pf.particles[p].y, gps_y));
    assert(areSame(pf.particles[p].theta, theta));
    assert(areSame(particle.weight, 1.));
    const int maxParticleId=lastParticleIdAtTestStart+numberOfParticles-numberOfParticlesToRemove;
    assert(particle.id<maxParticleId);// will the id's be preserved?
  }
  
  return 0;
}

int test() {
  testInit();
  testPredict();
  testError();
  testTransform();
  testDataAssociation();
  testParticleWeight();
  testFilterLandmarkMap();
  testResample();
  return 0;
}

int runFromFiles() {
  
  // parameters related to grading.
  int time_steps_before_lock_required = 100; // number of time steps before accuracy is checked by grader.
  double max_runtime = 45; // Max allowable runtime to pass [sec]
  double max_translation_error = 1; // Max allowable translation error to pass [m]
  double max_yaw_error = 0.05; // Max allowable yaw error [rad]
  
  // Start timer.
  int start = clock();
  
  //Set up parameters here
  double delta_t = 0.1; // Time elapsed between measurements [sec]
  double sensor_range = 50; // Sensor range [m]
  
  /*
   * Sigmas - just an estimate, usually comes from uncertainty of sensor, but
   * if you used fused data from multiple sensors, it's difficult to find
   * these uncertainties directly.
   */
  double sigma_pos [3] = {0.3, 0.3, 0.01}; // GPS measurement uncertainty [x [m], y [m], theta [rad]]
  double sigma_landmark [2] = {0.3, 0.3}; // Landmark measurement uncertainty [x [m], y [m]]
  
  // noise generation
  default_random_engine gen;
  normal_distribution<double> N_x_init(0, sigma_pos[0]);
  normal_distribution<double> N_y_init(0, sigma_pos[1]);
  normal_distribution<double> N_theta_init(0, sigma_pos[2]);
  normal_distribution<double> N_obs_x(0, sigma_landmark[0]);
  normal_distribution<double> N_obs_y(0, sigma_landmark[1]);
  double n_x, n_y, n_theta, n_range, n_heading;
  // Read map data
  Map map;
  if (!read_map_data("../runFromFiles/map_data.txt", map)) {
    cout << "Error: Could not open map file" << endl;
    return -1;
  }
  
  // Read position data
  vector<control_s> position_meas;
  if (!read_control_data("../runFromFiles/control_data.txt", position_meas)) {
    cout << "Error: Could not open position/control measurement file" << endl;
    return -1;
  }
  
  // Read ground truth data
  vector<ground_truth> gt;
  if (!read_gt_data("../runFromFiles/gt_data.txt", gt)) {
    cout << "Error: Could not open ground truth data file" << endl;
    return -1;
  }
  
  // Run particle filter!
  int num_time_steps = position_meas.size();
  ParticleFilter pf;
  double total_error[3] = {0,0,0};
  double cum_mean_error[3] = {0,0,0};
  
  for (int i = 0; i < num_time_steps; ++i) {
    cout << "Time step: " << i << endl;
    // Read in landmark observations for current time step.
    ostringstream file;
    file << "../runFromFiles/observation/observations_" << setfill('0') << setw(6) << i+1 << ".txt";
    vector<LandmarkObs> observations;
    if (!read_landmark_data(file.str(), observations)) {
      cout << "Error: Could not open observation file " << i+1 << endl;
      return -1;
    }
    
    // Initialize particle filter if this is the first time step.
    if (!pf.initialized()) {
      n_x = N_x_init(gen);
      n_y = N_y_init(gen);
      n_theta = N_theta_init(gen);
      pf.init(gt[i].x + n_x, gt[i].y + n_y, gt[i].theta + n_theta, sigma_pos);
    }
    else {
      // Predict the vehicle's next state (noiseless).
      pf.prediction(delta_t, sigma_pos, position_meas[i-1].velocity, position_meas[i-1].yawrate);
    }
    // simulate the addition of noise to noiseless observation data.
    vector<LandmarkObs> noisy_observations;
    LandmarkObs obs;
    for (int j = 0; j < observations.size(); ++j) {
      n_x = N_obs_x(gen);
      n_y = N_obs_y(gen);
      obs = observations[j];
      obs.x = obs.x + n_x;
      obs.y = obs.y + n_y;
      noisy_observations.push_back(obs);
    }
    
    // Update the weights and resample
    pf.updateWeights(sensor_range, sigma_landmark, noisy_observations, map);
    pf.resample();
    
    // Calculate and output the average weighted error of the particle filter over all time steps so far.
    vector<Particle> particles = pf.particles;
    int num_particles = particles.size();
    double highest_weight = 0.0;
    Particle best_particle;
    double weight_sum = 0.0;
    for (int i = 0; i < num_particles; ++i) {
      if (particles[i].weight > highest_weight) {
        highest_weight = particles[i].weight;
        best_particle = particles[i];
      }
      weight_sum += particles[i].weight;
    }
    
    cout << "highest w " << highest_weight << endl;
    cout << "average w " << weight_sum/num_particles << endl;
    
    double *avg_error = getError(gt[i].x, gt[i].y, gt[i].theta, best_particle.x, best_particle.y, best_particle.theta);
    
    for (int j = 0; j < 3; ++j) {
      total_error[j] += avg_error[j];
      cum_mean_error[j] = total_error[j] / (double)(i + 1);
    }
    
    // Print the cumulative weighted error
    cout << "Cumulative mean weighted error: x " << cum_mean_error[0] << " y " << cum_mean_error[1] << " yaw " << cum_mean_error[2] << endl;
    
    // If the error is too high, say so and then exit.
    if (i >= time_steps_before_lock_required) {
      if (cum_mean_error[0] > max_translation_error || cum_mean_error[1] > max_translation_error || cum_mean_error[2] > max_yaw_error) {
        if (cum_mean_error[0] > max_translation_error) {
          cout << "Your x error, " << cum_mean_error[0] << " is larger than the maximum allowable error, " << max_translation_error << endl;
        }
        else if (cum_mean_error[1] > max_translation_error) {
          cout << "Your y error, " << cum_mean_error[1] << " is larger than the maximum allowable error, " << max_translation_error << endl;
        }
        else {
          cout << "Your yaw error, " << cum_mean_error[2] << " is larger than the maximum allowable error, " << max_yaw_error << endl;
        }
        return -1;
      }
    }
  }
  
  // Output the runtime for the filter.
  int stop = clock();
  double runtime = (stop - start) / double(CLOCKS_PER_SEC);
  cout << "Runtime (sec): " << runtime << endl;
  
  // Print success if accuracy and runtime are sufficient (and this isn't just the starter code).
  if (runtime < max_runtime && pf.initialized()) {
    cout << "Success! Your particle filter passed!" << endl;
  }
  else if (!pf.initialized()) {
    cout << "This is the starter code. You haven't initialized your filter." << endl;
  }
  else {
    cout << "Your runtime " << runtime << " is larger than the maximum allowable runtime, " << max_runtime << endl;
    return -1;
  }
  
  return 0;
}

// for convenience
using json = nlohmann::json;

// Checks if the SocketIO event has JSON data.
// If there is data the JSON object in string format will be returned,
// else the empty string "" will be returned.
std::string hasData(std::string s) {
  auto found_null = s.find("null");
  auto b1 = s.find_first_of("[");
  auto b2 = s.find_first_of("]");
  if (found_null != std::string::npos) {
    return "";
  }
  else if (b1 != std::string::npos && b2 != std::string::npos) {
    return s.substr(b1, b2 - b1 + 1);
  }
  return "";
}

static bool RUNFROMFILES=false;
static bool RUNTESTS=false;

int main()
{
  if (RUNTESTS) {
    test();
  } else if (RUNFROMFILES) {
    return runFromFiles();
  } else {
    uWS::Hub h;
    
    //Set up parameters here
    double delta_t = 0.1; // Time elapsed between measurements [sec]
    double sensor_range = 50; // Sensor range [m]
    
    double sigma_pos [3] = {0.3, 0.3, 0.01}; // GPS measurement uncertainty [x [m], y [m], theta [rad]]
    double sigma_landmark [2] = {0.3, 0.3}; // Landmark measurement uncertainty [x [m], y [m]]
    
    // Read map data
    Map map;
    if (!read_map_data("../data/map_data.txt", map)) {
      cout << "Error: Could not open map file" << endl;
      return -1;
    }
    
    // Create particle filter
    ParticleFilter pf;
    
    h.onMessage([&pf,&map,&delta_t,&sensor_range,&sigma_pos,&sigma_landmark](uWS::WebSocket<uWS::SERVER> ws, char *data, size_t length, uWS::OpCode opCode) {
      // "42" at the start of the message means there's a websocket message event.
      // The 4 signifies a websocket message
      // The 2 signifies a websocket event
      
      std::cout << "data: <" << data << ">" << std::endl;
      if (length && length > 2 && data[0] == '4' && data[1] == '2')
      {
        
        auto s = hasData(std::string(data));
        if (s != "") {
          
          
          auto j = json::parse(s);
          std::string event = j[0].get<std::string>();
          
          if (event == "telemetry") {
            // j[1] is the data JSON object
            
            
            if (!pf.initialized()) {
              
              // Sense noisy position data from the simulator
              double sense_x = std::stod(j[1]["sense_x"].get<std::string>());
              double sense_y = std::stod(j[1]["sense_y"].get<std::string>());
              double sense_theta = std::stod(j[1]["sense_theta"].get<std::string>());
              
              pf.init(sense_x, sense_y, sense_theta, sigma_pos);
            }
            else {
              // Predict the vehicle's next state from previous (noiseless control) data.
              double previous_velocity = std::stod(j[1]["previous_velocity"].get<std::string>());
              double previous_yawrate = std::stod(j[1]["previous_yawrate"].get<std::string>());
              
              pf.prediction(delta_t, sigma_pos, previous_velocity, previous_yawrate);
            }
            
            // receive noisy observation data from the simulator
            // sense_observations in JSON format [{obs_x,obs_y},{obs_x,obs_y},...{obs_x,obs_y}]
            vector<LandmarkObs> noisy_observations;
            string sense_observations_x = j[1]["sense_observations_x"];
            string sense_observations_y = j[1]["sense_observations_y"];
            
            std::vector<float> x_sense;
            std::istringstream iss_x(sense_observations_x);
            
            std::copy(std::istream_iterator<float>(iss_x),
                      std::istream_iterator<float>(),
                      std::back_inserter(x_sense));
            
            std::vector<float> y_sense;
            std::istringstream iss_y(sense_observations_y);
            
            std::copy(std::istream_iterator<float>(iss_y),
                      std::istream_iterator<float>(),
                      std::back_inserter(y_sense));
            
            for(int i = 0; i < x_sense.size(); i++)
            {
              LandmarkObs obs;
              obs.x = x_sense[i];
              obs.y = y_sense[i];
              noisy_observations.push_back(obs);
            }
            
            // Update the weights and resample
            cout << "noisy_observations:" << noisy_observations.size() << std::endl;
            pf.updateWeights(sensor_range, sigma_landmark, noisy_observations, map);
            pf.resample();
            
            // Calculate and output the average weighted error of the particle filter over all time steps so far.
            vector<Particle> particles = pf.particles;
            int num_particles = particles.size();
            double highest_weight = -1.0;
            Particle best_particle;
            double weight_sum = 0.0;
            for (int i = 0; i < num_particles; ++i) {
              if (particles[i].weight > highest_weight) {
                highest_weight = particles[i].weight;
                best_particle = particles[i];
              }
              weight_sum += particles[i].weight;
            }
            cout << "highest w " << highest_weight << endl;
            cout << "average w " << weight_sum/num_particles << endl;
            
            json msgJson;
            msgJson["best_particle_x"] = best_particle.x;
            msgJson["best_particle_y"] = best_particle.y;
            msgJson["best_particle_theta"] = best_particle.theta;
            
            //Optional message data used for debugging particle's sensing and associations
            msgJson["best_particle_associations"] = pf.getAssociations(best_particle);
            msgJson["best_particle_sense_x"] = pf.getSenseX(best_particle);
            msgJson["best_particle_sense_y"] = pf.getSenseY(best_particle);
            
            auto msg = "42[\"best_particle\"," + msgJson.dump() + "]";
            // std::cout << msg << std::endl;
            ws.send(msg.data(), msg.length(), uWS::OpCode::TEXT);
            
          }
        } else {
          std::string msg = "42[\"manual\",{}]";
          ws.send(msg.data(), msg.length(), uWS::OpCode::TEXT);
        }
      }
      
    });
    
    // We don't need this since we're not using HTTP but if it's removed the program
    // doesn't compile :-(
    h.onHttpRequest([](uWS::HttpResponse *res, uWS::HttpRequest req, char *data, size_t, size_t) {
      const std::string s = "<h1>Hello world!</h1>";
      if (req.getUrl().valueLength == 1)
      {
        res->end(s.data(), s.length());
      }
      else
      {
        // i guess this should be done more gracefully?
        res->end(nullptr, 0);
      }
    });
    
    h.onConnection([&h](uWS::WebSocket<uWS::SERVER> ws, uWS::HttpRequest req) {
      std::cout << "Connected!!!" << std::endl;
    });
    
    h.onDisconnection([&h](uWS::WebSocket<uWS::SERVER> ws, int code, char *message, size_t length) {
      ws.close();
      std::cout << "Disconnected" << std::endl;
    });
    
    int port = 4567;
    if (h.listen(port))
    {
      std::cout << "Listening to port " << port << std::endl;
    }
    else
    {
      std::cerr << "Failed to listen to port" << std::endl;
      return -1;
    }
    h.run();
  }
}

