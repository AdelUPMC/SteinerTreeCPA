package algorithms;

import java.awt.Point;
import java.util.ArrayList;
import java.util.Collections;
import java.util.HashMap;
import java.util.List;
import java.util.Map;
import java.util.PriorityQueue;


public class DefaultTeam {

	public Tree2D calculSteiner(ArrayList<Point> points, int edgeThreshold, ArrayList<Point> hitPoints) {
		//étape 1: Construire le graphe pondéré complet  k.
		ArrayList<Edge> k = grapheK(points, hitPoints,edgeThreshold);

		//étape 2: Dans k, construire un arbre couvrant t0 de longueur totale des arêtes la plus petite possible.
		ArrayList<Edge> t0 = Kruskal(k,hitPoints);

		//étape 3: Dans t0, remplacer toute arête uv par un plus court chemin entre u et v dans G.
		ArrayList<Edge>h=replaceWithShortestPath(t0,points,edgeThreshold);

		//étape 4: Dans h, construire un arbre couvrant tprim de longueur totale des arêtes la plus petite possible.
		ArrayList<Edge> tprim =  Kruskal(h,points);

		//étape 5: Retourner l'arbre couvrant associé à tprim
		return edgesToTree(tprim,tprim.get(0).getP1());

	}
	public Tree2D calculSteinerBudget(ArrayList<Point> points, int edgeThreshold, ArrayList<Point> hitPoints) {
		double budget = 1664;
		Point initial = hitPoints.get(0);

		//étape 1: Construire le graphe pondéré complet  K 
		ArrayList<Edge> k = grapheK(points, hitPoints,edgeThreshold);

		//étape 2: Dans K, construire un arbre couvrant T0 de longeur totale des aretes la plus petite possible.
		ArrayList<Edge> t0 = Kruskal(k,hitPoints);

		//étape 3:récupérer la liste des arrêtes qui part du point initial(maison-mère) , couvrant le plus de hitPoints et satisfaisant la contrainte sur le budget.
		ArrayList<Edge> s= satisfyBudget( t0, initial, budget );

		//étape 4: Dans s, remplacer toute arête uv par un plus court chemin entre u et v dans G.
		ArrayList<Edge>budgetgraph=replaceWithShortestPath(s,points,edgeThreshold);

		// étape 5: Retourner l'arbre couvrant associé à budgetgraph
		return edgesToTree(budgetgraph,budgetgraph.get(0).getP1());
	}

	//------------------------------------------Utilitaires pour la version sans restriction budgétaire--------------------------------------------------------------------------------------------//


	public ArrayList<Edge> grapheK(ArrayList<Point> points, ArrayList<Point> hitpoints,int edgeThreshold) {
		//algorithme Floyd-warshall pour récupérer les matrices d'accès (plus courts chemins et distances)
		// int[][] path = calculShortestPaths( points, edgeThreshold );
		double[][] dist = calculDistances(points, edgeThreshold);
		double dij;

		ArrayList<Edge> k = new ArrayList<Edge>();

		//si i!=j on rajoute dans la liste des arêtes une arête de i vers j ayant pour poids la distance du plus court chemin entre i et j
		for (int i = 0;i<hitpoints.size();i++) {
			for (int j =0; j<hitpoints.size();j++) {
				if(i!=j) {
					dij= dist[points.indexOf(hitpoints.get(i))][points.indexOf(hitpoints.get(j))];
					k.add(new Edge(hitpoints.get(i),hitpoints.get(j),dij));
				}        
			}
		}
		return k;
	}

	public double[][] calculDistances(ArrayList<Point> points, int edgeThreshold) {
		double[][] dist =new double[points.size()][points.size()];
		double dij;
		for (int i = 0; i<dist.length; i++) {
			for (int j = 0;j<dist.length;j++) {
				//initialement la distance entre i et j correspond à la distance euclidienne entre les 2 points(si inférieure au seuil)
				dij= points.get(i).distance( points.get(j));
				if (dij <= edgeThreshold )dist[i][j] = dij;
				else dist[i][j] = Double.POSITIVE_INFINITY;
			}
		}

		for (int k = 0; k < dist.length; k++) {
			for (int i = 0; i < dist.length; i++ ) {
				for ( int j = 0; j < dist.length; j++ ) {
					if ( dist[i][j] > dist[i][k] + dist[k][j] ) {
						dist[i][j] = dist[i][k] + dist[k][j];
					}
				}
			}
		}
		return dist;
	}

	public ArrayList<Edge> Kruskal(ArrayList<Edge> edges, ArrayList<Point> vertices) {
		Collections.sort(edges);
		ArrayList<Edge> resultat = new ArrayList<>();
		Map<Point,Integer> etiquettes = initialiseLabels(vertices);		
		for ( Edge e : edges ) {
			//cas sans cycle: on ajoute l'arrête dans le résultat et on on met à jour les étiquettes
			if(!(etiquettes.get(e.getP1()) == etiquettes.get(e.getP2()))) {
				resultat.add(e);
				Integer labelp1 = (etiquettes.get(e.getP1()));
				Integer labelp2 = (etiquettes.get(e.getP2()));
				for (Edge edgeres : resultat) {
					if (etiquettes.get(edgeres.getP1()) == labelp1 ) {
						etiquettes.replace(edgeres.getP1(),labelp1, labelp2);
					}
					if (etiquettes.get(edgeres.getP2()) == labelp1) {
						etiquettes.replace(edgeres.getP2(),labelp1, labelp2);
					}
				}
			}     
		}
		return resultat;
	}



	public Map<Point, Integer> initialiseLabels(ArrayList<Point> points) {
		Map<Point, Integer> labels = new HashMap<Point, Integer>();
		for (int i=0;i<points.size();i++) {
			labels.put(points.get(i),i);
		}
		return labels;
	}
	public ArrayList<Edge> replaceWithShortestPath(ArrayList<Edge> t0,ArrayList<Point> points, int edgeThreshold){
		int[][] path = calculShortestPaths( points, edgeThreshold );
		ArrayList<Edge> resultat = new ArrayList<>();
		for (Edge e : t0) {
			int u = points.indexOf(e.getP1());
			int v = points.indexOf(e.getP2());
			int k = path[u][v];
			resultat.add(new Edge(e.getP1(),points.get(k)));
			while(k!=v){
				resultat.add(new Edge(points.get(k), points.get(path[k][v])));
				k= path[k][v];
			}
		}
		return resultat;

	}

	public int[][] calculShortestPaths(ArrayList<Point> points, int edgeThreshold) {
		int[][] paths=new int[points.size()][points.size()];
		for (int i=0;i<paths.length;i++) for (int j=0;j<paths.length;j++) paths[i][j]=i;

		double[][] dist=new double[points.size()][points.size()];

		for (int i=0;i<paths.length;i++) {
			for (int j=0;j<paths.length;j++) {
				if (i==j) {dist[i][i]=0; continue;}
				if (points.get(i).distance(points.get(j))<=edgeThreshold) dist[i][j]=points.get(i).distance(points.get(j)); else dist[i][j]=Double.POSITIVE_INFINITY;
				paths[i][j]=j;
			}
		}

		for (int k=0;k<paths.length;k++) {
			for (int i=0;i<paths.length;i++) {
				for (int j=0;j<paths.length;j++) {
					if (dist[i][j]>dist[i][k] + dist[k][j]){
						dist[i][j]=dist[i][k] + dist[k][j];
						paths[i][j]=paths[i][k];

					}
				}
			}
		}

		return paths;
	}	

	//------------------------------------------Utilitaires pour la version avec restriction budgétaire--------------------------------------------------------------------------------------------//
	/*
	 * implémentation de l'algorithme de Prim pour récupérer à chaque itération l'arrête de poids minimum
	 * en utilisant un tas minimum sous la forme d'une PriorityQueue<Edge>() avec un ordre sur le poids des arêtes.
	 */
	public ArrayList<Edge> satisfyBudget( ArrayList<Edge> edges, Point initial, double budget ) {
		ArrayList<Edge> resultat = new ArrayList<Edge>();
		Map<Point, PriorityQueue<Edge>> mapPointHeap = initialiseMapPointHeap(edges);
		//currentpoints: liste des points ajoutés déjà visités
		ArrayList<Point> currentpoints = new ArrayList<Point>();
		double used = 0.0;
		Edge edgeMin;
		currentpoints.add(initial);
		do {
			edgeMin = calculEdgeMin( mapPointHeap, currentpoints );
			if(used + edgeMin.getPoids() >= budget)break;
			mapPointHeap.get(edgeMin.getP1()).remove(edgeMin);
			mapPointHeap.get(edgeMin.getP2()).remove(edgeMin);
			if (!currentpoints.contains( edgeMin.getP1()))
				currentpoints.add( edgeMin.getP1());
			if (!currentpoints.contains( edgeMin.getP2()))
				currentpoints.add( edgeMin.getP2());
			used += edgeMin.getPoids();
			resultat.add( edgeMin);
		}while (used<budget);
		return resultat;
	}
	public Map<Point, PriorityQueue<Edge>> initialiseMapPointHeap( ArrayList<Edge> edges ) {
		Map<Point, PriorityQueue<Edge>> mapPointHeap = new HashMap<Point, PriorityQueue<Edge>>();
		for ( Edge edge : edges ) {
			//A chaque nouveau point visité, on initialise un tas min vide
			if ( !mapPointHeap.containsKey( edge.getP1() ) ) {
				mapPointHeap.put(edge.getP1(), new PriorityQueue<Edge>() );
			}
			if ( !mapPointHeap.containsKey( edge.getP2() ) ) {
				mapPointHeap.put(edge.getP2(), new PriorityQueue<Edge>() );
			}
			//si le point a déjà été visité on ajoute l'arête dans son tas min associé 
			mapPointHeap.get(edge.getP1()).add(edge);
			mapPointHeap.get(edge.getP2()).add(edge);

		}
		return mapPointHeap;
	}

	public Edge calculEdgeMin( Map<Point, PriorityQueue<Edge>> mapPointHeap, ArrayList<Point> currentpoints ) {
		/*On simule un tas min en utilisant une PriorityQueue<Edge>() avec un ordre sur le poids des arêtes, 
		 * mapPointHeap associe à chaque point le tas min qui contient toutes les aretes dans lesquelles 
		 * se trouve ce point.
		 * */
		PriorityQueue<Edge> candidates = new PriorityQueue<Edge>();
		for (Point p : currentpoints) {
			if(!(mapPointHeap.get(p).isEmpty()))
				candidates.add( mapPointHeap.get(p).peek());
		}
		//on retourne le premier élément du tas min pour avoir l'arête de poids minimum
		return candidates.peek();
	}

	//------------------------------------------Utilitaires---------------------------------------------------------------------------------------------------------------------------------------//
	public Tree2D edgesToTree(List<Edge> edges, Point root) {
		ArrayList<Edge> remainder = new ArrayList<Edge>();
		ArrayList<Point> subTreeRoots = new ArrayList<Point>();
		Edge current;
		while (edges.size()!=0) {
			current = edges.remove(0);
			if (current.getP1().equals(root)) {
				subTreeRoots.add(current.getP2());
			} else {
				if (current.getP2().equals(root)) {
					subTreeRoots.add(current.getP1());
				} else {
					remainder.add(current);
				}
			}
		}

		ArrayList<Tree2D> subTrees = new ArrayList<Tree2D>();
		for (Point subTreeRoot: subTreeRoots) subTrees.add(edgesToTree((ArrayList<Edge>)remainder.clone(),subTreeRoot));

		return new Tree2D(root, subTrees);
	}



}
