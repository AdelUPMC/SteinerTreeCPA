package algorithms;

import java.awt.Point;
public class Edge implements Comparable<Edge>{
		  private final Point p1;
		  private final Point p2;
		  private double poids;
		  
		  public Edge(Point p1, Point p2) {
		    this.p1 = p1;
		    this.p2 = p2;
		  }
		  public Edge(Point p1, Point p2, double poids) {
				this.p1 = p1;
				this.p2 = p2;
				this.poids = poids;
			}
		  public Point getP1() { return p1; }
		  public Point getP2() { return p2; }
		  
		  public double getPoids() {
			return poids;
		}
		@Override
			public int compareTo(Edge o) {
				int res = 0;
				if (this.getPoids() > o.getPoids())res = 1;
				else if (this.getPoids() < o.getPoids())res = -1;
				return res;
			}

		  @Override
		  public int hashCode() {
			  final int prime = 31;
			  int result = 1;
			  result = prime * result + ((p1 == null) ? 0 : p1.hashCode());
			  result = prime * result + ((p2 == null) ? 0 : p2.hashCode());
			  return result;
		  }

		  @Override
		  public boolean equals(Object obj) {
			  if (this == obj)
				  return true;
			  if (obj == null)
				  return false;
			  if (getClass() != obj.getClass())
				  return false;
			  Edge other = (Edge) obj;
			  if (p1 == null) {
				  if (other.p1 != null)
					  return false;
			  } else if (!p1.equals(other.p1))
				  return false;
			  if (p2 == null) {
				  if (other.p2 != null)
					  return false;
			  } else if (!p2.equals(other.p2))
				  return false;
			  return true;
		  }




}