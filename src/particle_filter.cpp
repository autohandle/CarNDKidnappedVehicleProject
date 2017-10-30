/*
 * particle_filter.cpp
 *
 *  Created on: Dec 12, 2016
 *      Author: Tiffany Huang
 */

#include <random>
#include <algorithm>
#include <iostream>
#include <numeric>
#include <math.h> 
#include <iostream>
#include <sstream>
#include <string>
#include <iterator>
#include <random>

#include "particle_filter.h"

using namespace std;

void ParticleFilter::init(const double x, const double y, const double theta, const double std[]) {
	// TODO: Set the number of particles. Initialize all particles to first position (based on estimates of 
	//   x, y, theta and their uncertainties from GPS) and all weights to 1. 
	// Add random Gaussian noise to each particle.
	// NOTE: Consult particle_filter.h for more information about this method (and others in this file).
  // This line creates a normal (Gaussian) distribution for x.]
  normal_distribution<double> dist_x(x, std[0]);
  
  // TODO: Create normal distributions for y and theta.
  normal_distribution<double> dist_y(y, std[1]);
  normal_distribution<double> dist_theta(theta, std[2]);
  
  for (int particle = 0; particle < num_particles; ++particle) {
    double sample_x, sample_y, sample_theta;
    
    // TODO: Sample  and from these normal distrubtions like this:
    // sample_x = dist_x(gen);
    // where "gen" is the random engine initialized earlier.
    sample_x=dist_x(gen);
    sample_y=dist_y(gen);
    sample_theta=dist_theta(gen);
    addParticleToFilter(createParticle(sample_x, sample_y, sample_theta));
  }

}

void ParticleFilter::prediction(const double delta_t, const double std_pos[]/* x,y,theta */, const double velocity, const double yaw_rate) {
	// TODO: Add measurements to each particle and add random Gaussian noise.
	// NOTE: When adding noise you may find std::normal_distribution and std::default_random_engine useful.
	//  http://en.cppreference.com/w/cpp/numeric/random/normal_distribution
	//  http://www.cplusplus.com/reference/random/default_random_engine/
  for (int p = 0; p < num_particles; ++p) {
    Particle& particle = particles[p];
    const double x0 = particle.x;
    const double y0 = particle.y;
    const double theta0 = particle.theta;
    const double xf = x0 +(velocity/yaw_rate)*(sin(theta0+yaw_rate*delta_t)-sin(theta0));
    const double yf = y0 +(velocity/yaw_rate)*(cos(theta0)-cos(theta0+yaw_rate*delta_t));
    const double thetaf = theta0+yaw_rate*delta_t;
    
    normal_distribution<double> noisyXf(xf, std_pos[0]);
    normal_distribution<double> noisyYf(yf, std_pos[1]);
    normal_distribution<double> noisyThetaf(thetaf, std_pos[2]);
    
    particle.x=noisyXf(gen);
    particle.y=noisyYf(gen);
    particle.theta=noisyThetaf(gen);
  }
}

void ParticleFilter::dataAssociation(std::vector<LandmarkObs>& theTransformedParticleObservations, std::vector<LandmarkObs>& theLandmarks) {
	// TODO: Find the predicted measurement that is closest to each observed measurement and assign the 
	//   observed measurement to this particular landmark.
	// NOTE: this method will NOT be called by the grading code. But you will probably find it useful to 
	//   implement this method and use it as a helper during the updateWeights phase.
  for (int p=0; p<theTransformedParticleObservations.size(); ++p) {
    LandmarkObs& prediction=theTransformedParticleObservations[p];
    prediction.id=theLandmarks[0].id;// assume 1st landmark is the best
    double bestDistance = dist(theLandmarks[0].x, theLandmarks[0].y, prediction.x, prediction.y);; //assume 1st landmark is the best
    for (int o=0; o<theLandmarks.size(); ++o) { // for every observation
      LandmarkObs& landmark = theLandmarks[o];
      const double distanceToLandmark=dist(landmark.x, landmark.y, prediction.x, prediction.y);
      if (distanceToLandmark<bestDistance) {
        // different observations can point to the same prediction
        prediction.id=landmark.id;// attach this prediction to this observation
        std::cout << "prediction:" << prediction.id << ", landmark:" << landmark.id << std::endl;
        bestDistance=distanceToLandmark;
      } // end if better distance
    }// end for all predicitions
  }// end for all observations
}

void ParticleFilter::updateWeights(double sensor_range, double std_landmark[], 
		const std::vector<LandmarkObs> &observations, const Map &map_landmarks) {
	// TODO: Update the weights of each particle using a mult-variate Gaussian distribution. You can read
	//   more about this distribution here: https://en.wikipedia.org/wiki/Multivariate_normal_distribution
	// NOTE: The observations are given in the VEHICLE'S coordinate system. Your particles are located
	//   according to the MAP'S coordinate system. You will need to transform between the two systems.
	//   Keep in mind that this transformation requires both rotation AND translation (but no scaling).
	//   The following is a good resource for the theory:
	//   https://www.willamette.edu/~gorr/classes/GeneralGraphics/Transforms/transforms2d.htm
	//   and the following is a good resource for the actual equation to implement (look at equation 
	//   3.33
	//   http://planning.cs.uiuc.edu/node99.html
}

void ParticleFilter::resample() {
	// TODO: Resample particles with replacement with probability proportional to their weight. 
	// NOTE: You may find std::discrete_distribution helpful here.
	//   http://en.cppreference.com/w/cpp/numeric/random/discrete_distribution

}

Particle ParticleFilter::SetAssociations(Particle particle, std::vector<int> associations, std::vector<double> sense_x, std::vector<double> sense_y)
{
	//particle: the particle to assign each listed association, and association's (x,y) world coordinates mapping to
	// associations: The landmark id that goes along with each listed association
	// sense_x: the associations x mapping already converted to world coordinates
	// sense_y: the associations y mapping already converted to world coordinates

	//Clear the previous associations
  // std::vector<int> associations;
  // std::vector<double> sense_x;
  // std::vector<double> sense_y;
	particle.associations.clear();
	particle.sense_x.clear();
	particle.sense_y.clear();

	particle.associations= associations;
 	particle.sense_x = sense_x;
 	particle.sense_y = sense_y;

 	return particle;
}

string ParticleFilter::getAssociations(Particle best)
{
	vector<int> v = best.associations;
	stringstream ss;
    copy( v.begin(), v.end(), ostream_iterator<int>(ss, " "));
    string s = ss.str();
    s = s.substr(0, s.length()-1);  // get rid of the trailing space
    return s;
}
string ParticleFilter::getSenseX(Particle best)
{
	vector<double> v = best.sense_x;
	stringstream ss;
    copy( v.begin(), v.end(), ostream_iterator<float>(ss, " "));
    string s = ss.str();
    s = s.substr(0, s.length()-1);  // get rid of the trailing space
    return s;
}
string ParticleFilter::getSenseY(Particle best)
{
	vector<double> v = best.sense_y;
	stringstream ss;
    copy( v.begin(), v.end(), ostream_iterator<float>(ss, " "));
    string s = ss.str();
    s = s.substr(0, s.length()-1);  // get rid of the trailing space
    return s;
}

const Particle& ParticleFilter::createParticle(const double theX, const double theY, const double theTheta) {
  Particle* const particle = new Particle;
  (*particle).id=particleId++;
  (*particle).x=theX;
  (*particle).y=theY;
  (*particle).theta=theTheta;
  return (*particle);
}

const void ParticleFilter::addParticleToFilter(const Particle& theParticle) {
  particles.push_back(theParticle);
}

const void ParticleFilter::deleteAllParticlesInFilter() {
  particles.clear();
}



