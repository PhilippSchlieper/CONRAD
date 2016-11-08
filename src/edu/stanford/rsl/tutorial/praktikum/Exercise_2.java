package edu.stanford.rsl.tutorial.praktikum;

import java.util.ArrayList;

import edu.stanford.rsl.conrad.data.numeric.Grid2D;
import edu.stanford.rsl.conrad.geometry.shapes.simple.Box;
import edu.stanford.rsl.conrad.geometry.shapes.simple.PointND;
import edu.stanford.rsl.conrad.geometry.shapes.simple.StraightLine;
import edu.stanford.rsl.conrad.geometry.transforms.Transform;
import edu.stanford.rsl.conrad.geometry.transforms.Translation;
import edu.stanford.rsl.conrad.numerics.SimpleOperators;
import edu.stanford.rsl.conrad.numerics.SimpleVector;
import ij.ImageJ;


public class Exercise_2 {
	
	public Grid2D computeSinogram(Grid2D grid, double angularRange, float detectorSize, float detectorSpacing, int numProj){
		double angularStep = angularRange / numProj;
		
		int maxDetectorIndex = (int) (detectorSize / detectorSpacing + 1);
		int maxAngularIndex = (int) (angularRange / angularStep + 1);
		
		final double samplingRate = 3.d; // # of samples per pixel
		Grid2D sino = new Grid2D(new float[maxAngularIndex*maxDetectorIndex], maxDetectorIndex, maxAngularIndex);
		sino.setSpacing(detectorSpacing, angularStep);

		// set up image bounding box in WC
		Translation trans = new Translation(
				-(grid.getSize()[0] * grid.getSpacing()[0])/2, -(grid.getSize()[1] * grid.getSpacing()[1])/2, -1);
		Transform inverse = trans.inverse();

		Box b = new Box((grid.getSize()[0] * grid.getSpacing()[0]), (grid.getSize()[1] * grid.getSpacing()[1]), 2);
		b.applyTransform(trans);

		for(int e=0; e<maxAngularIndex; ++e){
			// compute theta [rad] and angular functions.
			double theta = angularStep * e;
			double cosTheta = Math.cos(theta);
			double sinTheta = Math.sin(theta);

			for (int i = 0; i < maxDetectorIndex; ++i) {
				// compute s, the distance from the detector edge in WC [mm]
				double s = detectorSpacing * i - detectorSize / 2;
				// compute two points on the line through s and theta
				// We use PointND for Points in 3D space and SimpleVector for directions.
				PointND p1 = new PointND(s * cosTheta, s * sinTheta, .0d);
				PointND p2 = new PointND(-sinTheta + (s * cosTheta),
						(s * sinTheta) + cosTheta, .0d);
				// set up line equation
				StraightLine line = new StraightLine(p1, p2);
				// compute intersections between bounding box and intersection line.
				ArrayList<PointND> points = b.intersect(line);

				// only if we have intersections
				if (2 != points.size()){
					if(points.size() == 0) {
						line.getDirection().multiplyBy(-1.d);
						points = b.intersect(line);
					}
					if(points.size() == 0)
						continue;
				}

				PointND start = points.get(0); // [mm]
				PointND end = points.get(1);   // [mm]
				
				SimpleVector startVec = new SimpleVector(start.getAbstractVector());
				SimpleVector endVec = new SimpleVector(end.getAbstractVector());
				SimpleVector intVec = SimpleOperators.subtract(startVec, endVec);
				
				
			}
		}
	}
	
	
	
	public static void main(String[] args){
		
		ImageJ ij = new ImageJ();
		
		int phantomSize = 512;
		double angularRange = Math.PI;
		float detectorSize = 200;
		float detectorSpacing = 1.f;
		
		double [] phantSpacing = {detectorSpacing, detectorSpacing};
		Phantom phant = new Phantom(phantomSize, phantomSize, phantSpacing);
		
		phant.show();
		
	}

}
