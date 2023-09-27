package lineSimplification;

import java.awt.geom.Line2D;
import java.awt.geom.Point2D;
import java.io.File;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.HashSet;
import java.util.Iterator;
import java.util.List;
import java.util.Scanner;
import java.util.Set;
import java.util.StringTokenizer;

import gurobi.GRB;
import gurobi.GRBEnv;
import gurobi.GRBException;
import gurobi.GRBLinExpr;
import gurobi.GRBModel;
import gurobi.GRBVar;
import lineSimplification.DiGraph.DiGraphArc;
import lineSimplification.DiGraph.DiGraphNode;
import math.geom2d.Vector2D;

public class LineSimplification {

	/** 
	 * This method imports and reads a data from a file
	 * @param filepath - the file path of the data to be imported
	 * @param geographic - a boolean to determine whether the point data is geographic or not
	 * @return An ArrayList of Point 2D ArrayLists
	 * @throws Exception
	 */
	public static ArrayList<ArrayList<Point2D>> importFromFile(String filePath, boolean geographic) throws Exception{ 
		ArrayList<ArrayList<Point2D>> contours = new ArrayList<ArrayList<Point2D>>();
		Scanner s = new Scanner(new File(filePath));

		s.nextLine();

		while (s.hasNextLine()) {
			String first = s.nextLine();
			String[] t = first.split(";"); 

			StringTokenizer k = new StringTokenizer(t[0],"MULTILINESTRING, ()\"");

			ArrayList<Point2D> polyline = new ArrayList<Point2D>();
			while(k.hasMoreTokens()) {
				if (geographic) {
					Double x = Double.parseDouble(k.nextToken());
					Double y = Double.parseDouble(k.nextToken()); 
					Point2D p = LevelBasedProjection.WEBMERCATOR.fromLLtoPixel(new Point2D.Double(x, y), 5);
					polyline.add(p);
				} 
				else {
					Double x = Double.parseDouble(k.nextToken());
					Double y = Double.parseDouble(k.nextToken()); 
					Point2D p = new Point2D.Double(x, y);
					polyline.add(p);
				}
			}
			contours.add(polyline);

			StringTokenizer l = new StringTokenizer(t[2],",");
			Double elev = Double.parseDouble(l.nextToken());  //Get Elevation 
		}
		s.close();

		return contours;
	}



	/**
	 * This method implements the first task of Imai & Iri's line simplification algorithm 
	 * @param line - An ArrayList of Point 2D ArrayLists
	 * @param tolerance - Error tolerance
	 * @return A list of shortcut graphs
	 */
	public static ArrayList<DiGraph<Point2D, Line2D>> lineSimpl(ArrayList<ArrayList<Point2D>> line, double tolerance) {
		ArrayList<DiGraph<Point2D, Line2D>> listofShortcutGraph = new ArrayList<DiGraph<Point2D, Line2D>>() ;

		for (int p = 0; p < line.size(); p++) {
			ArrayList<Point2D> graph = line.get(p);
			DiGraph<Point2D,Line2D> shortcutGraph = new DiGraph<>();
			
			// Add the nodes of the shortcut graphs
			for (int i = 0; i < graph.size(); i++) {
				shortcutGraph.addNode(graph.get(i));
			}
			
			// Add the arc of the shortcut graphs based on the error tolerance
			for (int l = 0; l < shortcutGraph.n()-1; l++) {
				for (int j = l+1; j < shortcutGraph.n(); j++ ) {
					Line2D segment = new Line2D.Double(shortcutGraph.getNodeData(l), shortcutGraph.getNodeData(j));

					boolean add2Graph = true;
					for (int k = l + 1; k < j; k++ ) {
						double d = segment.ptLineDist(shortcutGraph.getNodeData(k));
						if (d > tolerance) {
							add2Graph = false;
						}
					}
					if (add2Graph) {
						shortcutGraph.addArc(shortcutGraph.getNode(l), shortcutGraph.getNode(j),segment);
					}	
				}
			}
			listofShortcutGraph.add(shortcutGraph);
		}
		return listofShortcutGraph;
	}

	
	
	/**
	 * This method implements a simple Integer Linear Program to simplify an ArrayList of polylines
	 * @param listofshortcutGraph - ArrayList of shortcut graphs
	 * @return An simplified graph
	 * @throws GRBException
	 * @throws NullPointerException
	 */
	public static ArrayList<DiGraph<Point2D, Line2D>> lineSimplOpt(ArrayList<DiGraph<Point2D, Line2D>> listofshortcutGraph) 
			throws GRBException,NullPointerException {
		ArrayList<DiGraph<Point2D, Line2D>> optShortcutGraph = new ArrayList<DiGraph<Point2D, Line2D>> ();

		ArrayList<DiGraphArc<Point2D, Line2D>> listOfAllArcs = new ArrayList<DiGraphArc<Point2D, Line2D>> ();
		
		ArrayList<Integer> graphIndex = new ArrayList<Integer> ();
		for (int u = 0; u < listofshortcutGraph.size(); u++) {
			DiGraph<Point2D, Line2D> graph = listofshortcutGraph.get(u);
			for (int z = 0; z < graph.m(); z++) {
				DiGraphArc<Point2D, Line2D> arcs = graph.getArc(z); //Get all the arcs in the graph
				listOfAllArcs.add(arcs);
				graphIndex.add(u);//Addition of the index of the graph that contains the arc 
			}
		}

		int l = listOfAllArcs.size(); //Number of arcs in the graph

		//Initialise the Gurobi model
		GRBModel m = new GRBModel(new GRBEnv());

		//Specify the variables to be solved
		List<GRBVar> vars = new ArrayList<GRBVar>();
	
		for (int f = 0; f < l; f++) {
			vars.add(m.addVar(0,1,0, GRB.BINARY, "x: ["+ f +"]")) ;
		}
		m.update();

		//Specify the objective function
		GRBLinExpr objFunx = new GRBLinExpr();
		for (int f = 0; f < l; f++ ) {
			objFunx.addTerm(1, vars.get(f));
		}
		m.setObjective(objFunx, GRB.MINIMIZE);


		//====== Constraints for solving the optimisation problem =======//
		for (int b = 0; b < listofshortcutGraph.size(); b++) {
			DiGraph<Point2D, Line2D> scgraph = listofshortcutGraph.get(b);

			///Constraint to ensure that only one arc leaves the source node 
			List<DiGraphArc<Point2D, Line2D>> outGoingArcsFromSource = scgraph.getNode(0).getOutgoingArcs();
			GRBLinExpr lh = new GRBLinExpr();
			for (int g = 0; g < outGoingArcsFromSource.size() ; g++) {
				int indexOfVar = getEdgeVarIndex(listOfAllArcs, outGoingArcsFromSource.get(g));
				lh.addTerm(1,vars.get(indexOfVar));
			}
			m.addConstr(lh, GRB.EQUAL, 1, "outgoingArcsSource = 1");

			///Constraint to ensure that only one arc gets to the target node 
			List<DiGraphArc<Point2D, Line2D>> inComingArcsToTarget = scgraph.getNode(scgraph.n() - 2).getIncomingArcs();
			GRBLinExpr lh1 = new GRBLinExpr();
			for (int h = 0; h < inComingArcsToTarget.size() ; h++) {
				int indexOfVar = getEdgeVarIndex(listOfAllArcs, inComingArcsToTarget.get(h));
				lh1.addTerm(1,vars.get(indexOfVar));
			}
			m.addConstr(lh1, GRB.EQUAL, 1, "incomingArcsTarget = 1");

			///Constraint to ensure that the number of arcs entering any intermediate node equals the number of arcs 
			///that leaves it
			for (int i = 1; i < scgraph.n()-1; i++) {
				List<DiGraphArc<Point2D, Line2D>> inComingArcs = scgraph.getNode(i).getIncomingArcs();
				List<DiGraphArc<Point2D, Line2D>> outGoingArcs = scgraph.getNode(i).getOutgoingArcs();

				GRBLinExpr linExpr = new GRBLinExpr(); 
				for (int j = 0; j < inComingArcs.size() ; j++) {
					int indexOfVar = getEdgeVarIndex(listOfAllArcs, inComingArcs.get(j));
					linExpr.addTerm(1, vars.get(indexOfVar));
				}

				for (int f = 0; f < outGoingArcs.size() ; f++) {
					int indexOfVar = getEdgeVarIndex(listOfAllArcs, outGoingArcs.get(f));
					linExpr.addTerm(-1, vars.get(indexOfVar));
				}
				m.addConstr(linExpr, GRB.EQUAL,0, "incomingArcs" + i + " = " + "outGoingArcs" + i);
			}
		}

		///Constraint to ensure that there are no intersections of arcs 
		for (int i = 0; i < l; i++) {
			for (int j = i+1; j < l; j++) {
				if (graphIndex.get(i) != graphIndex.get(j)) {
					if (listOfAllArcs.get(i).getArcData().intersectsLine(listOfAllArcs.get(j).getArcData())) { 
						GRBLinExpr interCon = new GRBLinExpr();
						interCon.addTerm(1.0, vars.get(i));
						interCon.addTerm(1.0, vars.get(j));
						m.addConstr(interCon, GRB.LESS_EQUAL, 1.0, "expr");
					}
				}
			}
		}

		m.optimize();
		//m.write("model.lp");


		//Printout of results
		for (int p = 0; p < l; p++) {
			double val = vars.get(p).get(GRB.DoubleAttr.X);
			DiGraph<Point2D,Line2D> optgraph = new DiGraph<>();
			if (Math.round(val) == 1) {
				optgraph.addArc(listOfAllArcs.get(p).getSource(), listOfAllArcs.get(p).getTarget(), listOfAllArcs.get(p).getArcData()); 
				optShortcutGraph.add(optgraph);
			}
		}
		return optShortcutGraph;
	}

	
	
	/**
	 * This method implements an advanced Integer Linear Program to simplify an ArrayList of polylines 
	 * @param groupofshortcutGraphs - ArrayList of shortcut graphs
	 * @param contours - ArrayList of Point 2D ArrayLists
	 * @param areaApprox - Fraction of the original area that final one should have
	 * @param thresholdAngle - Threshold angle
	 * @param fractionOfVertices - Fraction of the number of vertices
	 * @param edgeVarTradeOff - Trade_off coefficient for the edge variables
	 * @param angleVarTradeOff - Trade_off coefficient for the angle variables
	 * @param vertexVarTradeOff - Trade_off coefficient for the vertex variables
	 * @return An simplified graph
	 * @throws GRBException
	 * @throws NullPointerException
	 */
	public static ArrayList<DiGraph<Point2D, Line2D>> optimization(
			ArrayList<DiGraph<Point2D, Line2D>> groupofshortcutGraphs, 
			ArrayList<ArrayList<Point2D>> contours, 
			double areaApprox, 
			double thresholdAngle, 
			double fractionOfVertices,
			double edgeVarTradeOff, 
			double angleVarTradeOff,
			double vertexVarTradeOff) 
			throws GRBException,NullPointerException {
		
		ArrayList<DiGraph<Point2D, Line2D>> optimizedGraphs = new ArrayList<DiGraph<Point2D, Line2D>> ();
		
		ArrayList<DiGraph.DiGraphNode<Point2D, Line2D>> verticesList = new ArrayList<DiGraph.DiGraphNode<Point2D, Line2D>> ();

		ArrayList<DiGraph.DiGraphArc<Point2D, Line2D>> edgesList = new ArrayList<DiGraph.DiGraphArc<Point2D, Line2D>> ();
		ArrayList<Integer> graphIndexesofEdges = new ArrayList<Integer> ();
		
		ArrayList<Double> angleList = new ArrayList<Double>();
		ArrayList<Integer> indexesofAngles = new ArrayList<Integer>();
		
		for (int u = 0; u < groupofshortcutGraphs.size(); u++) {
			DiGraph<Point2D, Line2D> vertices = groupofshortcutGraphs.get(u);
			for (int z = 0; z < vertices.n(); z++) {
				DiGraph.DiGraphNode<Point2D, Line2D> vertex = vertices.getNode(z); 
				verticesList.add(vertex);
			}
		}
		
		for (int i = 0; i < groupofshortcutGraphs.size(); i++) {
			DiGraph<Point2D, Line2D> edges = groupofshortcutGraphs.get(i);
			for (int j = 0; j < edges.m(); j++) {
				DiGraph.DiGraphArc<Point2D, Line2D> edge = edges.getArc(j); 
				edgesList.add(edge);
				graphIndexesofEdges.add(i);
			}
		}
		
		int angleCount = 0;
		for (int i = 0; i < groupofshortcutGraphs.size(); i++) {
			DiGraph<Point2D, Line2D> simplifiedContour = groupofshortcutGraphs.get(i);
			for (int p = 1; p < simplifiedContour.n()-1; p++) {
				List<DiGraph.DiGraphArc<Point2D, Line2D>> incomingArcsOfNode = simplifiedContour.getNode(p).getIncomingArcs();
				List<DiGraph.DiGraphArc<Point2D, Line2D>> outgoingArcsOfNode = simplifiedContour.getNode(p).getOutgoingArcs();
				for (int j = 0; j < incomingArcsOfNode.size(); j++) {
					for (int k = 0; k < outgoingArcsOfNode.size(); k++) {
						double angle = calculateAngleDeg(incomingArcsOfNode.get(j), outgoingArcsOfNode.get(k));
						angleList.add(angle);
						int angleInd = angleCount++;
						indexesofAngles.add(angleInd);
					}
				}
			}
		}


		//The Gurobi model
		GRBModel grbModel = new GRBModel(new GRBEnv());

		//Variables to be solved
		/// Edge variables
		List<GRBVar> edgeVariables = new ArrayList<GRBVar>();
		for (int s = 0; s < edgesList.size(); s++) {
			edgeVariables.add(grbModel.addVar(0,1,0, GRB.BINARY, "x ("+ s +")")) ;
		}
		grbModel.update();
		
		///Angle variables
		List<GRBVar> angleVariables = new ArrayList<GRBVar>();
		for (int s = 0; s < angleList.size(); s++) {
			angleVariables.add(grbModel.addVar(0,1,0, GRB.BINARY, "y ("+ s +")")) ;
		}
		grbModel.update();
		
		/// Vertex variables
		List<GRBVar> vertexVariables = new ArrayList<GRBVar>();
		for (int s = 0; s < verticesList.size(); s++) {
			vertexVariables.add(grbModel.addVar(0,1,0, GRB.BINARY, "y ("+ s +")")) ;
		}
		grbModel.update();


		//======= Constraints for the optimisation problem ======//
		for (int i = 0; i < groupofshortcutGraphs.size(); i++) {
			DiGraph<Point2D, Line2D> ssimplifiedContour = groupofshortcutGraphs.get(i);
			
			///Constraint for the start vertex 
			List<DiGraph.DiGraphArc<Point2D, Line2D>> outGoingEdgesFromStVx = ssimplifiedContour.getNode(0).getOutgoingArcs();
			GRBLinExpr exp = new GRBLinExpr();
			for (int g = 0; g < outGoingEdgesFromStVx.size() ; g++) {
				int varIndex = getEdgeVarIndex(edgesList, outGoingEdgesFromStVx.get(g));
				exp.addTerm(1,edgeVariables.get(varIndex));
			}
			grbModel.addConstr(exp, GRB.EQUAL, 1, "outgoingEdgesFromStVx = 1");

			///Constraint for the last vertex 
			List<DiGraph.DiGraphArc<Point2D, Line2D>> inComingEdgesToLstVx = ssimplifiedContour.getNode(ssimplifiedContour.n() - 2).getIncomingArcs();
			GRBLinExpr exp1 = new GRBLinExpr();
			for (int h = 0; h < inComingEdgesToLstVx.size() ; h++) {
				int varIndex = getEdgeVarIndex(edgesList, inComingEdgesToLstVx.get(h));
				exp1.addTerm(1,edgeVariables.get(varIndex));
			}
			grbModel.addConstr(exp1, GRB.EQUAL, 1, "inComingEdgesToLstVx = 1");

			///Constraint to ensure that the number of edges entering any intermediate vertex equals the number of edges leaving it
			for (int k = 1; k < ssimplifiedContour.n()-1; k++) {
				List<DiGraph.DiGraphArc<Point2D, Line2D>> inComingArcs = ssimplifiedContour.getNode(k).getIncomingArcs();
				List<DiGraph.DiGraphArc<Point2D, Line2D>> outGoingArcs = ssimplifiedContour.getNode(k).getOutgoingArcs();

				GRBLinExpr linExpr = new GRBLinExpr();
				for (int j = 0; j < inComingArcs.size() ; j++) {
					int varIndex = getEdgeVarIndex(edgesList, inComingArcs.get(j));
					linExpr.addTerm(1, edgeVariables.get(varIndex));
				}

				for (int f = 0; f < outGoingArcs.size() ; f++) {
					int varIndex = getEdgeVarIndex(edgesList, outGoingArcs.get(f));
					linExpr.addTerm(-1, edgeVariables.get(varIndex));
				}
				grbModel.addConstr(linExpr, GRB.EQUAL,0, "incomingEdges" + i + " = " + "outGoingEdges" + i);
			}
			
			
			/// Vertex constraint
			GRBLinExpr expr = new GRBLinExpr();
			for (int j = 0; j < ssimplifiedContour.n()-1; j++) {
				int vertexId = getVertexVarIndex(verticesList, ssimplifiedContour.getNode(j));
				expr.addTerm(1, vertexVariables.get(vertexId));
			}
			grbModel.addConstr(expr, GRB.GREATER_EQUAL, fractionOfVertices*ssimplifiedContour.n(), null);
			
			
			/// Angle constraint - To check for acute angles
			for (int p = 1; p < ssimplifiedContour.n()-1; p++) {
				DiGraph.DiGraphNode<Point2D, Line2D> node = ssimplifiedContour.getNode(p);
				HashMap<Double,HashSet<DiGraph.DiGraphArc<Point2D, Line2D>>> setOfAnglesAndEdges = getAngle_ArcsList(node);

				Set<Double> Angles = setOfAnglesAndEdges.keySet();
				
				GRBLinExpr linExp = new GRBLinExpr();
				for (double key : Angles) {
					if(key <= thresholdAngle){
						double Angle = key;
						
						Iterator<DiGraph.DiGraphArc<Point2D, Line2D>> arcs = setOfAnglesAndEdges.get(key).iterator();
						int angleIndex = getAngleVarIndex(angleList, Angle);
						int incomingEdgeIndex = getEdgeVarIndex(edgesList, (DiGraph.DiGraphArc<Point2D, Line2D>) arcs.next());
						int outgoingEdgeIndex = getEdgeVarIndex(edgesList, (DiGraph.DiGraphArc<Point2D, Line2D>) arcs.next());
						linExp.addTerm(-1, angleVariables.get(angleIndex));
						linExp.addTerm(1, edgeVariables.get(incomingEdgeIndex));
						linExp.addTerm(1, edgeVariables.get(outgoingEdgeIndex));
					}
				}
				grbModel.addConstr(linExp, GRB.LESS_EQUAL,1, "incomingArcs" + i + " = " + "outGoingArcs" + i);	
			}
			
			
			/// Area constraint
			GRBLinExpr expArea = new GRBLinExpr();
			for (int k = 0; k < ssimplifiedContour.m()-1; k++) { 
				int startVertexId = getVertexVarIndex(verticesList, ssimplifiedContour.getArc(k).getSource());
				int endVertexId = getVertexVarIndex(verticesList, ssimplifiedContour.getArc(k).getTarget());
				
				if (endVertexId - startVertexId > 1) {
					ArrayList<Point2D> skippedNodes = new ArrayList<Point2D>();
					
					for (int y = startVertexId; y < endVertexId; y++) {
						skippedNodes.add(verticesList.get(y).getNodeData());
					}
					double area = areaofPolygon(skippedNodes);					
					expArea.addConstant(area);  
				}
			}
			grbModel.addConstr(expArea, GRB.GREATER_EQUAL, -areaApprox*areaofPolygon(contours.get(i)),"area constr"); 
		}	

		
		///No intersections of arcs constraint
		for (int m = 0; m < edgesList.size(); m++) { 
			for (int n = m+1; n < edgesList.size(); n++) { 
				if (graphIndexesofEdges.get(m) != graphIndexesofEdges.get(n)){
					if (edgesList.get(m).getArcData().intersectsLine(edgesList.get(n).getArcData())){ 
						GRBLinExpr consInte = new GRBLinExpr(); 
						consInte.addTerm(1, edgeVariables.get(m)); 
						consInte.addTerm(1, edgeVariables.get(n));
						grbModel.addConstr(consInte, GRB.LESS_EQUAL, 1, "ConstrOfInter");
					}
				}
			}
		}
		
		
		//The objective function
		GRBLinExpr obj = new GRBLinExpr();
		for (int t = 0; t < edgesList.size(); t++ ) {
			obj.addTerm(edgeVarTradeOff, edgeVariables.get(t));
		}
		
		for (int q = 0; q < angleList.size(); q++ ) {
			obj.addTerm(angleVarTradeOff, angleVariables.get(q));
		}
		
		for (int r = 0; r < verticesList.size(); r++ ) {
			obj.addTerm(vertexVarTradeOff, vertexVariables.get(r));
		}
		
		grbModel.setObjective(obj, GRB.MINIMIZE);

		grbModel.optimize();

		
		//Storing the results in an arrayList of graphs
		for (int pp = 0; pp < edgesList.size(); pp++) {
			double value = edgeVariables.get(pp).get(GRB.DoubleAttr.X);
			DiGraph<Point2D,Line2D> optimizedgraph = new DiGraph<>();
			if (Math.round(value) == 1) {
				optimizedgraph.addArc(edgesList.get(pp).getSource(), edgesList.get(pp).getTarget(), 
						edgesList.get(pp).getArcData()); 
				optimizedGraphs.add(optimizedgraph);
			}
		}
		return optimizedGraphs;
	}
	
	
	
	/**
	 * This method gets the index of a DiGraphArc in an ArrayList of DiGraphArcs
	 * @param arcs - An ArrayList of DiGraphArcs
	 * @param arc - A DiGraphArc
	 * @return An index
	 */
	public static int getEdgeVarIndex(ArrayList<DiGraphArc<Point2D, Line2D>> arcs, DiGraphArc<Point2D, Line2D> arc) {
		int index = 0;
		for (int i = 0; i < arcs.size(); i++) {
			if(arcs.get(i).equals(arc)) {
				index = i;
			}
		}
		return index;
	}
	
	
	
	/**
	 * This method gets the index of a DiGraphNode in an ArrayList of DiGraphNodes
	 * @param nodes - An ArrayList of DiGraphNodes
	 * @param node - A DiGraphNode
	 * @return An index
	 */
	public static int getVertexVarIndex(ArrayList<DiGraphNode<Point2D, Line2D>> nodes, DiGraphNode<Point2D, Line2D> node) {
		int index = 0;
		for (int i = 0; i < nodes.size(); i++) {
			if(nodes.get(i).equals(node)) {
				index = i;
			}
		}
		return index;
	}


	
	/**
	 * This method gets the index of an angle in an ArrayList of angles
	 * @param angles - An ArrayList of angles
	 * @param angle - An angle
	 * @return An index
	 */
	public static int getAngleVarIndex(ArrayList<Double> angles, Double angle) {
		int ind = 0;
		for (int i = 0; i < angles.size(); i++) {
			if(angles.get(i).equals(angle)) {
				ind = i;
			}
		}
		return ind;
	}
	
	
	
	/**
	 * This method computes the area  of a polygon using the coordinates of the vertices
	 * @param contour - An ArrayList of 2D points
	 * @return The area of the polygon
	 */
	public static double areaofPolygon(ArrayList<Point2D> contour){
		double area = 0.0;
		for (int i = 0; i < contour.size()-1; i++) { 
			area += (contour.get(i).getX() * contour.get(i+1).getY()) - (contour.get(i+1).getX() * contour.get(i).getY());
		}
		return Math.abs(area / 2.0);
	}
	
	

	/**
	 * This method computes the area of a closed graph using the coordinates of the nodes
	 * @param graph -  A closed graph
	 * @return The area of the graph
	 */
	public static double areaOfClosedGraph(DiGraph<Point2D, Line2D> graph){
		double area = 0.0;
        for (int i = 0; i < graph.n()-1; i++) 
            area += (graph.getNodeData(i).getX() * graph.getNodeData(i+1).getY()) - (graph.getNodeData(i+1).getX() * graph.getNodeData(i).getY());
        return Math.abs(area / 2.0);
	}
	
	
	
	/**
	 * This method computes the angle between two vectors in degrees
	 * @param incomingArc - An incoming vector
	 * @param outgoingArc - An outgoing vector
	 * @return An angle
	 */
	public static double calculateAngleDeg(DiGraphArc<Point2D, Line2D> incomingArc, DiGraphArc<Point2D, Line2D> outgoingArc) {
		double incX1 = incomingArc.getSource().getNodeData().getX();
		double incX2 = incomingArc.getTarget().getNodeData().getX();
		double incY1 = incomingArc.getSource().getNodeData().getY();
		double incY2 = incomingArc.getTarget().getNodeData().getY();
		
		double outX2 = outgoingArc.getTarget().getNodeData().getX();
		double outX1 = outgoingArc.getSource().getNodeData().getX();
		double outY2 = outgoingArc.getTarget().getNodeData().getY();
		double outY1 = outgoingArc.getSource().getNodeData().getY();
		Vector2D vectorOfIncomingEdge = new Vector2D(incX1 - incX2, incY1 - incY2);
		Vector2D vectorOfOutgoingEdge = new Vector2D(outX2 - outX1, outY2 - outY1);

		double angle = Math.acos(Vector2D.dot(vectorOfIncomingEdge,vectorOfOutgoingEdge) /
				(vectorOfIncomingEdge.norm() * vectorOfOutgoingEdge.norm())) * 180/Math.PI;
		
		return angle;
	}
	
	
	
	/**
	 * This method gets a HashMap of an angle and a HashSet of the DiGraphArcs forming the angle
	 * @param vertex - A DiGraphNode
	 * @return A HashMap of an angle and a HashSet of arcs
	 */
	public static HashMap<Double,HashSet<DiGraph.DiGraphArc<Point2D, Line2D>>> getAngle_ArcsList(
			DiGraph.DiGraphNode<Point2D, Line2D> vertex){
		HashMap<Double,HashSet<DiGraph.DiGraphArc<Point2D, Line2D>>> listOfAngles = 
				new HashMap<Double,HashSet<DiGraph.DiGraphArc<Point2D, Line2D>>>();

		List<DiGraph.DiGraphArc<Point2D, Line2D>> incomingArcsOfNode = vertex.getIncomingArcs();
		List<DiGraph.DiGraphArc<Point2D, Line2D>> outgoingArcsOfNode = vertex.getOutgoingArcs();
		HashSet<DiGraph.DiGraphArc<Point2D, Line2D>> arcsFormingAngle = new HashSet<DiGraph.DiGraphArc<Point2D, Line2D>>();
		
		for (int j = 0; j < incomingArcsOfNode.size(); j++) {
			for (int k = 0; k < outgoingArcsOfNode.size(); k++) {
				double angle = calculateAngleDeg(incomingArcsOfNode.get(j), outgoingArcsOfNode.get(k));
				arcsFormingAngle.add(incomingArcsOfNode.get(j));
				arcsFormingAngle.add(outgoingArcsOfNode.get(k));
				listOfAngles.put(angle, arcsFormingAngle);
			}
		}
		return listOfAngles;
	}

}




