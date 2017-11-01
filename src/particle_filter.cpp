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
  
  if (num_particles==0) {
    num_particles=500;
  }
  assert(num_particles>0);
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

  is_initialized=true;
}

bool isYawRateZero(const double theValue) {
  return isZero(theValue);
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
    
    double xf;
    double yf;
    double thetaf;
    
    if (isYawRateZero(yaw_rate)) {
      const double deltaV = delta_t*velocity;
      xf = x0 + deltaV*cos(theta0);
      yf = y0 + deltaV*sin(theta0);
      thetaf = theta0;
      //cout << "x0:" << x0 << ", xf:" << xf << ", velocity:" << velocity << ", deltaV:" << (velocity*cos(theta0)) << std::endl;
      //cout << "y0:" << y0 << ", yf:" << yf << ", theta0:" << theta0 << ", deltaV:" << (velocity*sin(theta0)) << std::endl;
    } else {
      xf = x0 +(velocity/yaw_rate)*(sin(theta0+yaw_rate*delta_t)-sin(theta0));
      yf = y0 +(velocity/yaw_rate)*(cos(theta0)-cos(theta0+yaw_rate*delta_t));
      thetaf = theta0+yaw_rate*delta_t;
    }
    
    normal_distribution<double> noisyXf(xf, std_pos[0]);
    normal_distribution<double> noisyYf(yf, std_pos[1]);
    normal_distribution<double> noisyThetaf(thetaf, std_pos[2]);
    //normal_distribution<double> noisyXf(xf, 0.);
    //normal_distribution<double> noisyYf(yf, 0.);
    //normal_distribution<double> noisyThetaf(thetaf, 0.);
    
    particle.x=noisyXf(gen);
    particle.y=noisyYf(gen);
    assert(!std::isnan(particle.y) && !std::isnan(particle.y));
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
    prediction.id=0;// assume 1st landmark is the best
    double bestDistance = dist(theLandmarks[0].x, theLandmarks[0].y, prediction.x, prediction.y);; //assume 1st landmark is the best
    for (int lm=0; lm<theLandmarks.size(); ++lm) { // for every observation
      LandmarkObs& landmark = theLandmarks[lm];
      const double distanceToLandmark=dist(landmark.x, landmark.y, prediction.x, prediction.y);
      if (distanceToLandmark<bestDistance) {
        // different prediction can point to the same landmark
        prediction.id=lm;// attach this prediction to this landmark
        //std::cout << "prediction:" << prediction.id << ", landmark:" << landmark.id << std::endl;
        bestDistance=distanceToLandmark;
      } // end if better distance
    }// end for all predicitions
  }// end for all observations
}

bool isParicleWeightZero(const double theWeight) {
  return isZero(theWeight, 1.e-100);
}

bool isParicleWeightGreaterThanZero(const Particle& theParticle) {
  return !isParicleWeightZero(theParticle.weight);
}

const double totalParticleWeight(std::vector<Particle>& theParticles) {
  double totalWeight=0.;
  for (int p=0; p<theParticles.size();p++) {
    Particle& particle = theParticles[p];
    if (isParicleWeightGreaterThanZero(particle)) {
      totalWeight+=particle.weight;
    }
  }
  return totalWeight;
}

const double totalParticleWeight(ParticleFilter& theParticleFilter) {
  return totalParticleWeight(theParticleFilter.particles);
}

void ParticleFilter::updateWeights(double sensor_range, double std_landmark[], 
		const std::vector<LandmarkObs>& observations, const Map& map_landmarks) {
  //cout << "observations:" << observations.size() << std::endl;

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
  assert(particles.size()>0);
  for (int p=0; p<particles.size();p++) {
    Particle& particle = particles[p];
    std::vector<LandmarkObs> /* xMap, yMap */  observationsInMapCoordinates=ParticleFilter::transform(particle, observations);

    std::vector<LandmarkObs> landmarksInSensorRange=filterLandmarkMap(map_landmarks, particle, sensor_range);
    dataAssociation(observationsInMapCoordinates, landmarksInSensorRange); // link observations to landmarks
    const double particleWeight = ParticleFilter::particleWeight(observationsInMapCoordinates, landmarksInSensorRange, std_landmark);// if 1 landmark is way off (e.g. exp-infinity) then particleWeight will be 0
    //cout << "observations:" << observations.size() << ", observationsInMapCoordinates:" << observationsInMapCoordinates.size() << ", weight:" << particle.weight << "->" << particleWeight << std::endl;
    particle.weight=particleWeight;
    //cout << "particle["+to_string(p)+"].weight: "+to_string(particle.weight) << std::endl;
  }
  //cout << "observations:" << observations.size() << std::endl;
  
  if (!(totalParticleWeight(particles)>0.)) {
    Particle& particle = particles[0];
    std::vector<LandmarkObs> /* xMap, yMap */  observationsInMapCoordinates=ParticleFilter::transform(particle, observations);
    
    std::vector<LandmarkObs> landmarksInSensorRange=filterLandmarkMap(map_landmarks, particle, sensor_range);
    dataAssociation(observationsInMapCoordinates, landmarksInSensorRange); // link observations to landmarks
    const double particleWeight = ParticleFilter::particleWeight(observationsInMapCoordinates, landmarksInSensorRange, std_landmark);
  }
  
  assert(totalParticleWeight(particles)>0.);
}

void ParticleFilter::resample() {
	// TODO: Resample particles with replacement with probability proportional to their weight. 
	// NOTE: You may find std::discrete_distribution helpful here.
	//   http://en.cppreference.com/w/cpp/numeric/random/discrete_distribution
  std::vector<Particle> particlesToResample;
  // first collect particles with a nonzero weight
  double totalWeight=0.;
  for (int p=0; p<particles.size();p++) {
    Particle& particle = particles[p];
    if (isParicleWeightGreaterThanZero(particle)) {
      totalWeight+=particle.weight;
      particlesToResample.push_back(particle);
    }
  }
  
  assert(totalWeight>0.);
  assert(particlesToResample.size()>0);
  //const int particleShortage=num_particles-particlesToResample.size();
  //cout << "particleShortage:" << particleShortage << std::endl;
  /*
  // then add some random particles based on the ones woth nonzero weight
  for (int p=0; p<particleShortage; p++) {
    Particle& particle = particlesToResample[p%particlesToResample.size()];
    const Particle& randomClone = createRandomParticle(particle);
    particlesToResample.push_back(randomClone);
    totalWeight+=randomClone.weight;
  }
  
  assert(particlesToResample.size()==num_particles);
  */
  //double cummulativeParticleProbability[particlesWithNonZeroWeight.size()];
  std::vector<double> cummulativeParticleProbability;
  double cummulativeProbability=0.;
  for (int p=0; p<particlesToResample.size();p++) {
    const Particle& particle = particlesToResample[p];
    const double particleProbabilty=particle.weight/totalWeight;
    cummulativeProbability+=particleProbabilty;
    //cummulativeParticleProbability[p]=cummulativeProbability;
    cummulativeParticleProbability.push_back(cummulativeProbability);

  }
  assert(areSame(cummulativeProbability, 1.0));
  assert(cummulativeParticleProbability.size()>0);
  std::uniform_real_distribution<double> distribution(0.,1.);
  std::vector<Particle> resamples;
  for (int p=0; p<num_particles;p++) {// forall particle slots in the original number of samples
    const double randomProbability = distribution(gen);
    int cpi=0;
    while (cummulativeParticleProbability[cpi]<randomProbability) {
      cpi+=1;
      //cout << "randomProbability:" << randomProbability << ", cummulativeParticleProbability[" << cpi << "]:" << cummulativeParticleProbability[cpi] << std::endl;
    }
    resamples.push_back(particlesToResample[cpi]);
  }
  particles=resamples;
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

int& ParticleFilter::currentParticleId() {
  return particleId;
}

const int ParticleFilter::nextParticleId() {
  return currentParticleId()++;
}

const Particle& ParticleFilter::createParticle(const double theX, const double theY, const double theTheta) {
  Particle* const particle = new Particle;
  //(*particle).id=particleId++;
  (*particle).id=nextParticleId();
  (*particle).x=theX;
  (*particle).y=theY;
  (*particle).theta=theTheta;
  (*particle).weight=1.0;// all particles are equally likely
  return (*particle);
}

const Particle& ParticleFilter::createRandomParticle(const double theX, const double theY, const double theTheta, const double theWeight) {
  const double std[] = {0.3, 0.3, 0.1};
  
  Particle* const particle = new Particle;
  
  normal_distribution<double> randomX(theX, std[0]);
  normal_distribution<double> randomY(theY, std[1]);
  normal_distribution<double> randomTheta(theTheta, std[2]);

  //(*particle).id=particleId++;
  (*particle).id=-nextParticleId();// flag as a randomw particle
  (*particle).x=randomX(gen);
  (*particle).y=randomY(gen);
  (*particle).theta=randomTheta(gen);
  assert(!isZero(theWeight, 1.e-100));
  (*particle).weight=theWeight;
  return (*particle);
}

const Particle& ParticleFilter::createRandomParticle(const Particle& theParticle) {
  return createRandomParticle(theParticle.x, theParticle.y, theParticle.theta, theParticle.weight/2.);
}


const void ParticleFilter::addParticleToFilter(const Particle& theParticle) {
  particles.push_back(theParticle);
}

const void ParticleFilter::deleteAllParticlesInFilter() {
  particles.clear();
}



