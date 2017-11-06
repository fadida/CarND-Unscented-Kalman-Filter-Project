# Unscented Kalman Filter Project
Self-Driving Car Engineer Nanodegree Program

This project works with Udacity self driving simulator (for term 2) and
uses unscented Kalman filter in order to predict vehicle position from
radar and lidar sensors.

## The algorithm
The first step in the algorithm is taking the initial measurements in order to initialize the Kalman filter.

After initialization, the alogrithm is built to predict first and then update the location and state using the current measurements using a different update method for each sensor.

The predict method creates a matrix of sigma points which is used by the update method of each sensor in order to compare the measurement with the filter prediction. 

After the measurements update the program calculates the RMSE in order to analyze the accuracy of the algorithm.

Also added, is the option to calculate the print the filter consistency using the NIS test for each sensor.
This test was added in order to help tweaking the filter prediction noise marix.


