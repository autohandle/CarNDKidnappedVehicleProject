/*
 * particle_filter.h
 *
 * 2D particle filter class.
 *  Created on: Dec 12, 2016
 *      Author: Tiffany Huang
 */

#ifndef PARTICLE_FILTER_H_
#define PARTICLE_FILTER_H_

#include "helper_functions.h"
#include <random>
#include "assert.h"

struct Particle {

	int id;
	double x;
	double y;
	double theta;
	double weight;
	std::vector<int> associations;
	std::vector<double> sense_x;
	std::vector<double> sense_y;
};



class ParticleFilter {
	
	// Number of particles to draw
	int num_particles; 
	
	
	
	// Flag, if filter is initialized
	bool is_initialized;
	
	// Vector of weights of all particles
	std::vector<double> weights;
	
public:
  
  static LandmarkObs& /* xMap, yMap */ transform(const double particleX, const double particleY, const double theRoationfromCarToParticle, const double observationX, const double observationY) {
    LandmarkObs* const inMapCoordinates = new LandmarkObs;
    (*inMapCoordinates).id=-1;// flag landmark as unlinked
    (*inMapCoordinates).x = particleX+cos(theRoationfromCarToParticle)*observationX-sin(theRoationfromCarToParticle)*observationY;
    (*inMapCoordinates).y = particleY+sin(theRoationfromCarToParticle)*observationX+ cos(theRoationfromCarToParticle)*observationY;
    return (*inMapCoordinates);
  }
  
  static LandmarkObs& /* xMap, yMap */ transform(const Particle& theParticle, const LandmarkObs& theObservation) {
    return transform(theParticle.x, theParticle.y, theParticle.theta, theObservation.x, theObservation.y);
  }
	
  static std::vector<LandmarkObs> /* xMap, yMap */ transform(const Particle theParticle, std::vector<LandmarkObs> theObservations) {
    std::vector<LandmarkObs> transformedObservations;
    for (int o=0; o<theObservations.size(); o++) {
      const LandmarkObs& observation = theObservations[o];
      const LandmarkObs& transformedObservation = transform(theParticle, observation);
      assert(!std::isnan(transformedObservation.x) && !std::isnan(transformedObservation.y));
      transformedObservations.push_back(transformedObservation);
    }
    return transformedObservations;
  }
  
  static const double particleWeight(const LandmarkObs& theLandmark/* u */, const LandmarkObs& theTransformedObservation /* x,y */, const double theObservationMeasurementSigmas[]/* sx, sy */) {
    const double sigmaX=theObservationMeasurementSigmas[0];
    const double sigmaXSquared = sigmaX*sigmaX;
    const double sigmaY=theObservationMeasurementSigmas[1];
    const double sigmaYSquared = sigmaY*sigmaY;
    const double xDiff=theTransformedObservation.x-theLandmark.x;
    const double xDiffSquared=xDiff*xDiff;
    const double yDiff=theTransformedObservation.y-theLandmark.y;
    const double yDiffSquared=yDiff*yDiff;
    
    return (1./(2*M_PI*sigmaX*sigmaY))*exp(-1.*((xDiffSquared/(2.*sigmaXSquared)) + (yDiffSquared/(2.*sigmaYSquared))) );
  }

  static const double particleWeight(std::vector<LandmarkObs>& theTransformedObservations /* x,y */, std::vector<LandmarkObs>& theLandmarks/* u */, const double theObservationMeasurementSigmas[]/* sx, sy */) {
    //assert(theTransformedObservations.size()<=theLandmarks.size()); data has more observations than landmarks
    double cummulativeWeight=1.0;
    for (int obs=0; obs<theTransformedObservations.size(); obs++) {
      const LandmarkObs& observation = theTransformedObservations[obs];
      assert((observation.id>=0) && (observation.id<theLandmarks.size()));// observation.id==-1 -> not linked to a landmark
      const LandmarkObs& landmark = theLandmarks[observation.id];
      const double landmarkObservationWeight = particleWeight(landmark, observation, theObservationMeasurementSigmas);
      cummulativeWeight*=landmarkObservationWeight;// if 1 landmark is way off (e.g. exp-infinity) then cummulativeWeight will be 0
    }
    return cummulativeWeight;
  }
  
  static std::vector<LandmarkObs> filterLandmarkMap(const Map& map_landmarks, const Particle& theParticle, const double theSensorRange) {
    std::vector<LandmarkObs> landmarksInSensorRange;
    for (int lm=0; lm<map_landmarks.landmark_list.size(); lm++) {
      const Map::single_landmark_s& landmarkInMap = map_landmarks.landmark_list[lm];
      if (fabs(fabs(landmarkInMap.x_f)-fabs(theParticle.x))<theSensorRange && fabs(fabs(landmarkInMap.y_f)-fabs(theParticle.y))<theSensorRange) {
        LandmarkObs landmarkInRange;
        landmarkInRange.id=landmarkInMap.id_i;
        landmarkInRange.x=landmarkInMap.x_f;
        landmarkInRange.y=landmarkInMap.y_f;
        landmarksInSensorRange.push_back(landmarkInRange);
      }
    }
    assert(landmarksInSensorRange.size()>0);
    return landmarksInSensorRange;
  }
  
  // Set of current particles
	std::vector<Particle> particles;

	// Constructor
	// @param num_particles Number of particles
	ParticleFilter() : num_particles(0), is_initialized(false), gen(rd()) {}

	// Destructor
	~ParticleFilter() {}

	/**
	 * init Initializes particle filter by initializing particles to Gaussian
	 *   distribution around first position and all the weights to 1.
	 * @param x Initial x position [m] (simulated estimate from GPS)
	 * @param y Initial y position [m]
	 * @param theta Initial orientation [rad]
	 * @param std[] Array of dimension 3 [standard deviation of x [m], standard deviation of y [m]
	 *   standard deviation of yaw [rad]]
	 */
	void init(const double x, const double y, const double theta, const double std[]);

	/**
	 * prediction Predicts the state for the next time step
	 *   using the process model.
	 * @param delta_t Time between time step t and t+1 in measurements [s]
	 * @param std_pos[] Array of dimension 3 [standard deviation of x [m], standard deviation of y [m]
	 *   standard deviation of yaw [rad]]
	 * @param velocity Velocity of car from t to t+1 [m/s]
	 * @param yaw_rate Yaw rate of car from t to t+1 [rad/s]
	 */
	void prediction(const double delta_t, const double std_pos[], const double velocity, const double yaw_rate);
	
	/**
	 * dataAssociation Finds which observations correspond to which landmarks (likely by using
	 *   a nearest-neighbors data association).
	 * @param predicted Vector of predicted landmark observations
	 * @param observations Vector of landmark observations
	 */
	void dataAssociation(std::vector<LandmarkObs>& predicted, std::vector<LandmarkObs>& observations);
	
	/**
	 * updateWeights Updates the weights for each particle based on the likelihood of the 
	 *   observed measurements. 
	 * @param sensor_range Range [m] of sensor
	 * @param std_landmark[] Array of dimension 2 [Landmark measurement uncertainty [x [m], y [m]]]
	 * @param observations Vector of landmark observations
	 * @param map Map class containing map landmarks
	 */
	void updateWeights(double sensor_range, double std_landmark[], const std::vector<LandmarkObs> &observations,
			const Map &map_landmarks);
	
	/**
	 * resample Resamples from the updated set of particles to form
	 *   the new set of particles.
	 */
	void resample();

	/*
	 * Set a particles list of associations, along with the associations calculated world x,y coordinates
	 * This can be a very useful debugging tool to make sure transformations are correct and assocations correctly connected
	 */
	Particle SetAssociations(Particle particle, std::vector<int> associations, std::vector<double> sense_x, std::vector<double> sense_y);
	
	std::string getAssociations(Particle best);
	std::string getSenseX(Particle best);
	std::string getSenseY(Particle best);

	/**
	 * initialized Returns whether particle filter is initialized yet or not.
	 */
	const bool initialized() const {
		return is_initialized;
	}
  
  const void setNumberOfParticles(const int theNumberOfParticles) {
    if (initialized()) {
      throw std::logic_error("Already Initialized");
    }
    num_particles=theNumberOfParticles;
  }

  int& currentParticleId();
  const int nextParticleId();
  const Particle& createParticle(const double theX, const double theY, const double theTheta);
  const Particle& createRandomParticle(const double theX, const double theY, const double theTheta, const double theWeight);
  const Particle& createRandomParticle(const Particle& theParticle);
  const void addParticleToFilter(const Particle& theParticle);
  const void deleteAllParticlesInFilter();
 
private:
  
  //std::default_random_engine gen;
  //std::default_random_engine generator;
  std::random_device rd;
  std::mt19937 gen;

};

static int particleId=0;

#endif /* PARTICLE_FILTER_H_ */
