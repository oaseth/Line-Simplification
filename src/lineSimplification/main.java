package lineSimplification;

import java.awt.Color;
import java.awt.geom.Line2D;
import java.awt.geom.Path2D;
import java.awt.geom.Point2D;
import java.io.IOException;
import java.util.ArrayList;

import javax.swing.JFrame;

import gurobi.GRBException;


public class main {

	public static void main(String[] args) throws Exception{
		long timebegun = System.currentTimeMillis(); //Get the start time
		ArrayList<ArrayList<Point2D>> contours = LineSimplification.importFromFile("Contours1.csv", false);
		
		//Creation of the simplified contours
		double errorTolerance = 5;
		double thresholdAngle = 100; // Angle in degrees
		double fractionOfVertices = 1;
		double areaApprox = 0.01;
		double edgeVarTradeOff = 0.3; 
		double angleVarTradeOff = 0.2;  
		double vertexVarTradeOff = 0.3; 
				
		ArrayList<DiGraph<Point2D, Line2D>> simplifiedGraph = LineSimplification.optimization(
				LineSimplification.lineSimpl(contours,errorTolerance), 
				contours, 
				thresholdAngle,
				fractionOfVertices,
				areaApprox, 
				edgeVarTradeOff, 
				angleVarTradeOff,
				vertexVarTradeOff);
		
		//Visualisation of the original contours and the simplified ones
		Viewer view = new Viewer();
		JFrame frame = new JFrame("LINE SIMPLIFICATION    [Original contour - BLACK; Simplified contour - MAGENTA] "
				+ "TOLERANCE = " + errorTolerance + ", MINIMUM ANGLE = " + thresholdAngle);
		frame.setSize(1050,700);
		frame.add(view);
		
		///Loop through all the polylines of original contours and add them to the Viewer object, "view"
		for (int i = 0; i < contours.size(); i++) {
			ArrayList<Point2D> contour = contours.get(i);
			Path2D ct = new Path2D.Double();
	
			ct.moveTo(contour.get(0).getX(), contour.get(0).getY());//Starting point of the path

			for (int j = 1; j < contour.size(); j++) {
				ct.lineTo(contour.get(j).getX(), contour.get(j).getY());//Addition of the remaining points to the path
				view.addShape(ct,Color.BLACK,1.5);	
			}
		}
		
		///Loop through all the graphs of simplified contours and add them to the Viewer object, "view"
		for (int s = 0; s < simplifiedGraph.size(); s++) {
			DiGraph<Point2D, Line2D> simplifiedContour = simplifiedGraph.get(s);
			for (int j = 0; j < simplifiedContour.m(); j++) {
				Line2D arc = simplifiedContour.getArc(j).getArcData();
				view.addShape(arc,Color.MAGENTA,1.5);
			}
		}
		
		view.moveToCenter();
		frame.setVisible(true);
		frame.setDefaultCloseOperation(JFrame.EXIT_ON_CLOSE);
		
		long timestopped = System.currentTimeMillis(); //Get the stop time
		long timespan = (timestopped - timebegun);    //Get the running time the algorithm
		System.out.println("Running time: "+ timespan + " milliseconds");  //Print
	}
}
