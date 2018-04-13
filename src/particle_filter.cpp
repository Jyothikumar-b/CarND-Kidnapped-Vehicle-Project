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

#include "particle_filter.h"

using namespace std;

void ParticleFilter::init(double x, double y, double theta, double std[]) {
	default_random_engine gen;
	double std_x, std_y, std_theta; // Standard deviations for x, y, and theta
	std_x=std[0];
	std_y=std[1];
	std_theta=std[2];
	
	// TODO: Set the number of particles.
	num_particles=100;
	
	//Initialize all particles to first position (based on estimates of 
	//   x, y, theta and their uncertainties from GPS) and all weights to 1.
	normal_distribution<double> dist_x(x, std_x);
	normal_distribution<double> dist_y(y, std_y);
	normal_distribution<double> dist_theta(theta, std_theta);
	
	for (int i = 0; i < num_particles; ++i) {
		Particle p;
		
		p.x = dist_x(gen);
		p.y = dist_y(gen);
		p.theta = dist_theta(gen);
		p.weight=1.0;
		weights.push_back(1.0);
		particles.push_back(p);
		
	}
	
	is_initialized = true;
	
}

void ParticleFilter::prediction(double delta_t, double std_pos[], double velocity, double yaw_rate) {
	// TODO: Add measurements to each particle and add random Gaussian noise.
	// NOTE: When adding noise you may find std::normal_distribution and std::default_random_engine useful.
	//  http://en.cppreference.com/w/cpp/numeric/random/normal_distribution
	//  http://www.cplusplus.com/reference/random/default_random_engine/
	
	default_random_engine gen;
	
	normal_distribution<double> dist_x(0, std_pos[0]);
	normal_distribution<double> dist_y(0, std_pos[1]);
	normal_distribution<double> dist_theta(0, std_pos[2]);
	
	for(int i=0;i<num_particles;++i){
		
		if( fabs(yaw_rate) < 0.0001){
            //if velocity is too close to zero
            particles[i].x += velocity * delta_t * cos(particles[i].theta) + dist_x(gen);
            particles[i].y += velocity * delta_t * sin(particles[i].theta) + dist_y(gen);
        }else{
			particles[i].x = particles[i].x + (velocity/yaw_rate) * ( sin(particles[i].theta + yaw_rate*delta_t) - sin(particles[i].theta) ) + dist_x(gen);
			particles[i].y = particles[i].y + (velocity/yaw_rate) * ( cos(particles[i].theta) - cos(particles[i].theta + yaw_rate*delta_t) ) + dist_y(gen);
			particles[i].theta = particles[i].theta + yaw_rate*delta_t + dist_theta(gen);
		}
	}

}

void ParticleFilter::dataAssociation(std::vector<LandmarkObs> predicted, std::vector<LandmarkObs>& observations) {
	// TODO: Find the predicted measurement that is closest to each observed measurement and assign the 
	//   observed measurement to this particular landmark.
	// NOTE: this method will NOT be called by the grading code. But you will probably find it useful to 
	//   implement this method and use it as a helper during the updateWeights phase.
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
	
	for(int i_p=0;i_p<num_particles;i_p++){
		
		std::vector<LandmarkObs> SelectedLandMarks;
		std::vector<double> Landmark_dist;
		
		//Landmark within the given range and save the distance (check each particle)
		for(int i_L=0;i_L<int(map_landmarks.landmark_list.size());i_L++){
			double distance=dist(map_landmarks.landmark_list[i_L].x_f,map_landmarks.landmark_list[i_L].y_f,particles[i_p].x,particles[i_p].y);
			if (distance<sensor_range){
				LandmarkObs selectedLM;
				selectedLM.x=map_landmarks.landmark_list[i_L].x_f;
				selectedLM.y=map_landmarks.landmark_list[i_L].y_f;
				selectedLM.id=map_landmarks.landmark_list[i_L].id_i;
				SelectedLandMarks.push_back(selectedLM);
				Landmark_dist.push_back(distance);
			}
		}
		
		//1. convert observation value into global map for that particle
		std::vector<LandmarkObs> ObservedLandMarks;
		for(int obs=0;obs<int(observations.size());obs++){
			LandmarkObs detectedLM;
			detectedLM.x = particles[i_p].x + observations[obs].x*cos(particles[i_p].theta) - observations[obs].y*sin(particles[i_p].theta);
			detectedLM.y = particles[i_p].y + observations[obs].x*sin(particles[i_p].theta) + observations[obs].y*cos(particles[i_p].theta);
			ObservedLandMarks.push_back(detectedLM);
		}
		
		
		//2. find the best pair with probability/weight
		particles[i_p].weight=1.0;
		for(int obs=0;obs<int(ObservedLandMarks.size());obs++){
			double distance=0,cdist=0;
			int min=-1;
			for(int act=0;act<int(SelectedLandMarks.size());act++){
				cdist=dist(SelectedLandMarks[act].x,SelectedLandMarks[act].y,ObservedLandMarks[obs].x,ObservedLandMarks[obs].y);
				if(cdist<distance || act==0){
					distance=cdist;
					min=act;
				}
			}
			double x = SelectedLandMarks[min].x - ObservedLandMarks[obs].x;
			double y = SelectedLandMarks[min].y - ObservedLandMarks[obs].y;

			//3. Calculate total probability of the particle
			particles[i_p].weight=particles[i_p].weight*( 1/(2*M_PI*std_landmark[0]*std_landmark[1]) * exp(-0.5*((x*x)/(std_landmark[0]*std_landmark[0])+
																												(y*y)/(std_landmark[1]*std_landmark[1])))); 
			
		}
		weights[i_p]=particles[i_p].weight;
		
	}
	
	
}

void ParticleFilter::resample() {
	// TODO: Resample particles with replacement with probability proportional to their weight. 
	// NOTE: You may find std::discrete_distribution helpful here.
	//   http://en.cppreference.com/w/cpp/numeric/random/discrete_distribution
	default_random_engine gen;

	vector<Particle> resampled;

	discrete_distribution<int> dist_w(weights.begin(), weights.end());

	for (auto &particle : particles){

		resampled.push_back(particles[dist_w(gen)]);

	}

	particles = resampled;
	
}

Particle ParticleFilter::SetAssociations(Particle& particle, const std::vector<int>& associations, 
                                     const std::vector<double>& sense_x, const std::vector<double>& sense_y)
{
    //particle: the particle to assign each listed association, and association's (x,y) world coordinates mapping to
    // associations: The landmark id that goes along with each listed association
    // sense_x: the associations x mapping already converted to world coordinates
    // sense_y: the associations y mapping already converted to world coordinates

    particle.associations= associations;
    particle.sense_x = sense_x;
    particle.sense_y = sense_y;
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
