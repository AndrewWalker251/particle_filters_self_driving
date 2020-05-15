/**
 * particle_filter.cpp
 *
 * Created on: Dec 12, 2016
 * Author: Tiffany Huang
 */

#include "particle_filter.h"

#include <math.h>
#include <algorithm>
#include <iostream>
#include <iterator>
#include <numeric>
#include <random>
#include <string>
#include <vector>


#include "helper_functions.h"

using std::string;
using std::vector;
using std::normal_distribution;
using std::discrete_distribution;

void ParticleFilter::init(double x, double y, double theta, double std[]) {
  /**
   * TODO: Set the number of particles. Initialize all particles to 
   *   first position (based on estimates of x, y, theta and their uncertainties
   *   from GPS) and all weights to 1. 
   * TODO: Add random Gaussian noise to each particle.
   * NOTE: Consult particle_filter.h for more information about this method 
   *   (and others in this file).
   */
  num_particles = 50;  // TODO: Set the number of particles
  
  // this will be used to randomly select from the distribution.
  std::default_random_engine gen;
  
  // we need to randomly create the y,x and the theta from separate distributions
  // use the measement as the mean.
  normal_distribution<double> dist_x(x, std[0]);
  normal_distribution<double> dist_y(y, std[1]);
  normal_distribution<double> dist_theta(theta, std[2]);

  for  (int N = 0; N <num_particles; N++)
  {
   Particle particle;
   particle.id = N;
   particle.x = dist_x(gen);
   particle.y = dist_y(gen);
   particle.theta = dist_theta(gen);
   particle.weight = 1.0;

   particles.push_back(particle);
   weights.push_back(1.0);
  }
  is_initialized = true;
}

void ParticleFilter::prediction(double delta_t, double std_pos[], 
                                double velocity, double yaw_rate) {
  /**
   * TODO: Add measurements to each particle and add random Gaussian noise.
   * NOTE: When adding noise you may find std::normal_distribution 
   *   and std::default_random_engine useful.
   *  http://en.cppreference.com/w/cpp/numeric/random/normal_distribution
   *  http://www.cplusplus.com/reference/random/default_random_engine/
   */
  std::default_random_engine gen2;
  // for each particle move the designated amount. 
  for (int i =0; i < num_particles; i++)
  {
    float x = particles[i].x;
    float y = particles[i].y;
    float theta = particles[i].theta;
    float yrdt = yaw_rate*delta_t;
    
  // Need to check for division by 0
    
   if (fabs(yaw_rate) < 0.00001)
   {
    x = x + velocity*delta_t *cos(theta);
    y = y + velocity*delta_t * sin(theta);
    theta = theta;
   }
   else
   {
    x = x + (velocity/yaw_rate)*(sin(theta + yrdt) - sin(theta));
    y = y + (velocity/yaw_rate)*(cos(theta) - cos(theta + yrdt));
    theta = theta + yrdt;
   }  
     
    normal_distribution<double> dist_x(x, std_pos[0]);
    normal_distribution<double> dist_y(y, std_pos[1]);
    normal_distribution<double> dist_theta(theta, std_pos[2]);
    
    particles[i].x = dist_x(gen2);
    particles[i].y = dist_y(gen2);
    particles[i].theta = dist_theta(gen2);
  }
}

void ParticleFilter::dataAssociation(vector<LandmarkObs> predicted, 
                                     vector<LandmarkObs>& observations) {
  /**
   * TODO: Find the predicted measurement that is closest to each 
   *   observed measurement and assign the observed measurement to this 
   *   particular landmark.
   * NOTE: this method will NOT be called by the grading code. But you will 
   *   probably find it useful to implement this method and use it as a helper 
   *   during the updateWeights phase.
   */
  
}
void ParticleFilter::updateWeights(double sensor_range, double std_landmark[], 
                                   const vector<LandmarkObs> &observations, 
                                   const Map &map_landmarks) {
  /**
   * TODO: Update the weights of each particle using a mult-variate Gaussian 
   *   distribution. You can read more about this distribution here: 
   *   https://en.wikipedia.org/wiki/Multivariate_normal_distribution
   * NOTE: The observations are given in the VEHICLE'S coordinate system. 
   *   Your particles are located according to the MAP'S coordinate system. 
   *   You will need to transform between the two systems. Keep in mind that
   *   this transformation requires both rotation AND translation (but no scaling).
   *   The following is a good resource for the theory:
   *   https://www.willamette.edu/~gorr/classes/GeneralGraphics/Transforms/transforms2d.htm
   *   and the following is a good resource for the actual equation to implement
   *   (look at equation 3.33) http://planning.cs.uiuc.edu/node99.html
   */
  // Using guidance from provided Q&A.
  
  //Step 1. Convert observations from car coordinates to map coordinates at the particle.
  //vector<double> weights;
  for (int i=0; i<num_particles ; i++)
  	{
    	vector<double> temporary_weights_holder;
        vector<int> associations;
	    vector<double> sense_x;
	    vector<double> sense_y;
  		float x_part = particles[i].x;
  		float y_part = particles[i].y;
  		float theta = particles[i].theta;

 		// I'm assuming that the car always points upwards so the theta we use to do the tranpose is theta.
  
 		// For each of the sensors
    	for (unsigned  O=0; O<observations.size() ; O++)
        {	
          	float x_obs = observations[O].x;
  			float y_obs = observations[O].y;	
          	// convert observation to the paricle location and position.
  			double x_map;
  			x_map = x_part + (cos(theta) * x_obs) - (sin(theta) * y_obs);
  			double y_map;
  			y_map = y_part + (sin(theta) * x_obs) + (cos(theta) * y_obs);
  
  			// Now we're in the correct coordinates which real object is the closest
  
  			// for all of the real objects
  			// keep track of which one is closest 
            // what if none are in range - thenwe have a weight 0. 
            int closest = -1;
            double closest_distance = 10000000;
          	for (unsigned  r = 0; r< map_landmarks.landmark_list.size(); r++)
            {
               double distance = dist(x_map, y_map, map_landmarks.landmark_list[r].x_f, map_landmarks.landmark_list[r].y_f);
               //std::cout << "distance" << distance << std::endl;
               //std::cout << "closest" << closest_distance << std::endl;
               if (distance < closest_distance)
               {
                 closest_distance = distance;
                 closest = r;
            	}
            }
          
			// for each observation take the closest lankmark and calculate p.
              
            if (closest != -1)
            {
                float nearest_real_x = map_landmarks.landmark_list[closest].x_f;
            	float nearest_real_y = map_landmarks.landmark_list[closest].y_f;
                int landmark_id = map_landmarks.landmark_list[closest].id_i;
          
                // Calculate the weight.   
                // calculate normalization term
                double gauss_norm;
                gauss_norm = 1 / (2 * M_PI * std_landmark[0] * std_landmark[1]);

                // calculate exponent
                double exponent;
                exponent = (pow((x_map - nearest_real_x), 2) / (2 * pow(std_landmark[0], 2)))
                             + (pow((y_map - nearest_real_y), 2) / (2 * pow(std_landmark[1], 2)));

                // calculate weight using normalization terms and exponent
                double weight;
                weight = gauss_norm * exp(-exponent);
                //std::cout << "y_map " << y_map << std::endl;
                //std::cout << "nearest_y " << nearest_real_y << std::endl;
                //std::cout << "x_map " << x_map << std::endl;
                //std::cout << "nearest_x " << nearest_real_x << std::endl;
				//std::cout << "exponent " << exponent << std::endl;
                //std::cout << "weight " << weight << std::endl;
                // add this weight to the list for this particle. Once you're got them all multiply together. 
                temporary_weights_holder.push_back(weight);
              	associations.push_back(landmark_id);
				sense_x.push_back(x_map);
				sense_y.push_back(y_map);
                
            }
          //else
          //{
          //  temporary_weights_holder.push_back(0);
          //}
        } 
        // Once you're all the weights for the particle multiply them together.
		  SetAssociations(particles[i], associations, sense_x, sense_y);
              
		  //std::cout << "particle id " << i << std::endl;
    	  for(unsigned  out=0; out<temporary_weights_holder.size(); ++out)
          {
  			//std::cout << temporary_weights_holder[out] << ' ';
          }
          float final_weight = temporary_weights_holder[0];

          for(unsigned  v = 1; v< temporary_weights_holder.size(); v++)
          {
          final_weight = final_weight * temporary_weights_holder[v];
          }
        
    
    particles[i].weight = final_weight;       
    //std::cout << "final weight " << final_weight << std::endl;
    //std::cout << "particle weight " << particles[i].weight << std::endl;
    weights[i] = (particles[i].weight);
    } 
    std::cout <<  "weights array ";
    for(unsigned  out=0; out<weights.size(); ++out)
          {
  			//std::cout << weights[out] << ' ';
          }

    
  
  
}     
void ParticleFilter::resample() {
  /**
   * TODO: Resample particles with replacement with probability proportional 
   *   to their weight. 
   * NOTE: You may find std::discrete_distribution helpful here.
   *   http://en.cppreference.com/w/cpp/numeric/random/discrete_distribution
   */
  std::default_random_engine gen3;
  std::cout <<  "weights array 2 ";
    for(unsigned  out=0; out<weights.size(); ++out)
          {
  			//std::cout << weights[out] << ' ';
          }
  discrete_distribution<int> distribution(weights.begin(),weights.end());
  
  vector<Particle> resample_particles; 
  
  for (int i=0; i<num_particles; i++)
  {
    resample_particles.push_back(particles[distribution(gen3)]);
  }
 std::cout << std::endl << "particles before";
  for(unsigned  out=0; out<particles.size(); ++out)
          {
  			//std::cout << particles[out].id << ' ';
          }
  
  particles = resample_particles;
  
  std::cout << std::endl << "particles selected ";
  for(unsigned  out=0; out<particles.size(); ++out)
          {
  			//std::cout << particles[out].id << ' ';
          }
  
  
}

void ParticleFilter::SetAssociations(Particle& particle, 
                                     const vector<int>& associations, 
                                     const vector<double>& sense_x, 
                                     const vector<double>& sense_y) {
  // particle: the particle to which assign each listed association, 
  //   and association's (x,y) world coordinates mapping
  // associations: The landmark id that goes along with each listed association
  // sense_x: the associations x mapping already converted to world coordinates
  // sense_y: the associations y mapping already converted to world coordinates
  particle.associations= associations;
  particle.sense_x = sense_x;
  particle.sense_y = sense_y;
}

string ParticleFilter::getAssociations(Particle best) {
  vector<int> v = best.associations;
  std::stringstream ss;
  copy(v.begin(), v.end(), std::ostream_iterator<int>(ss, " "));
  string s = ss.str();
  s = s.substr(0, s.length()-1);  // get rid of the trailing space
  return s;
}

string ParticleFilter::getSenseCoord(Particle best, string coord) {
  vector<double> v;

  if (coord == "X") {
    v = best.sense_x;
  } else {
    v = best.sense_y;
  }

  std::stringstream ss;
  copy(v.begin(), v.end(), std::ostream_iterator<float>(ss, " "));
  string s = ss.str();
  s = s.substr(0, s.length()-1);  // get rid of the trailing space
  return s;
}