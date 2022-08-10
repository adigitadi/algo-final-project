package algo;

import java.awt.geom.Point2D;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collections;
import java.util.HashMap;
import java.util.HashSet;
import java.util.List;
import java.util.Map;
import java.util.PriorityQueue;
import java.util.Set;

class Vertex {
  String name;
  double x;
  double y;

  public Vertex(String name, double x, double y) {
    this.name = name;
    this.x = x;
    this.y = y;
  }

  @Override
  public String toString() {
    return name;
  }

  @Override
  public boolean equals(Object o) {
    if(o == this) {
      return true;
    }
    if(!(o instanceof Vertex)) {
      return false;
    }
    Vertex v = (Vertex) o;

    return Double.compare(x, v.x) == 0 && Double.compare(y, v.y)== 0 && name.equalsIgnoreCase(v.name); 
  }

}

class Edge {
  Vertex source;
  Vertex destination;
  double distance;

  public Edge(Vertex source, Vertex destination, double distance) {
    this.source = source;
    this.destination = destination;
    this.distance = distance;
  }
}

public class Kmeans{

  public static void main(String []args){
    System.out.println("Hello, World!");

    Vertex v2 = new Vertex("a", 1, 1);
    Vertex v3 = new Vertex("b", 2, 3);
    Vertex v4 = new Vertex("c", -2, 2);
    Vertex v5 = new Vertex("d", -1, 4);

    Vertex v6 = new Vertex("e", -1, -1);
    Vertex v7 = new Vertex("f", 2, -3);
    Vertex v8 = new Vertex("g", -2, -1);
    Vertex v9 = new Vertex("h", 5, 0);

    Vertex v10 = new Vertex("i", 1, 1);
    Vertex v11 = new Vertex("j", 2, 3);
    Vertex v12 = new Vertex("k", -2, 2);
    Vertex v13 = new Vertex("l", -1, 4);

    Vertex v14 = new Vertex("m", 1, 1);
    Vertex v15 = new Vertex("n", 2, 3);
    Vertex v16 = new Vertex("o", -2, 2);
    Vertex v17 = new Vertex("p", -1, 4);

    List<Vertex> vertices = new ArrayList<>(Arrays.asList(v2,v3,v4,v5,v6,v7,v8,v9,v10,v11,v12,v13,v14,v15,v16,v17));

    int clusters = 4;

    int perCar = 5;

    Map<String, List<Vertex>> group = cluster(vertices, clusters, createCentroids(clusters, vertices), perCar);

    Map<String, List<Vertex>> group1 = cluster(vertices, clusters, reiterate(group, vertices), perCar);

    while(!group1.equals(group)) {
      group = group1;
      group1 = cluster(vertices, clusters, reiterate(group, vertices), perCar);
      System.out.println(group.entrySet());
    }

    Vertex s = new Vertex("s",0, 0);

    List<List<String>> allPaths = new ArrayList<>();

    findAllPaths(allPaths, group1, s);


  }

  private static void findAllPaths(List<List<String>> allPaths, Map<String, List<Vertex>> group1, Vertex startingPoint) {
    List<Integer> tour = new ArrayList<>();
    double minTourCost = Double.POSITIVE_INFINITY;
    for(List<Vertex> list : group1.values()) {
      double[][] distance = new double[list.size() + 1][list.size() + 1];
      list.add(0, startingPoint);

      for(int i =0; i < list.size() -1; i++) {
        for(int j =i; j < list.size(); j++) {
          distance[i][j] = getDistance(list.get(i).x, list.get(i).y, list.get(j).x, list.get(j).y);
          distance[j][i] = distance[i][j];
        }
      }

      //start
      int n = distance.length;
      int finalState = (int) Math.pow(2, n) - 1;
      double[][] memo = new double[n][(int) Math.pow(2, n)];

      int start = 0;
      // Add all outgoing edges from the starting node to memo table.
      for (int end = 0; end < n; end++) {
        if (end == start) continue;
        memo[end][(1 << start) | (1 << end)] = distance[start][end];
      }

      for (int r = 3; r <= n; r++) {
        for (int subset : combinations(r, n)) {
          if (notIn(start, subset)) continue;
          for (int next = 0; next < n; next++) {
            if (next == start || notIn(next, subset)) continue;
            int subsetWithoutNext = subset ^ (1 << next);
            double minDist = Double.POSITIVE_INFINITY;
            for (int end = 0; end < n; end++) {
              if (end == start || end == next || notIn(end, subset)) continue;
              double newDistance = memo[end][subsetWithoutNext] + distance[end][next];
              if (newDistance < minDist) {
                minDist = newDistance;
              }
            }
            memo[next][subset] = minDist;
          }
        }
      }

      // Connect tour back to starting node and minimize cost.
      for (int i = 0; i < n; i++) {
        if (i == start) continue;
        double tourCost = memo[i][finalState] + distance[i][start];
        if (tourCost < minTourCost) {
          minTourCost = tourCost;
        }
      }

      int lastIndex = start;
      int state = finalState;
      tour.add(start);

      // Reconstruct TSP path from memo table.
      for (int i = 1; i < n; i++) {

        int index = -1;
        for (int j = 0; j < n; j++) {
          if (j == start || notIn(j, state)) continue;
          if (index == -1) index = j;
          double prevDist = memo[index][state] + distance[index][lastIndex];
          double newDist  = memo[j][state] + distance[j][lastIndex];
          if (newDist < prevDist) {
            index = j;
          }
        }

        tour.add(index);
        state = state ^ (1 << index);
        lastIndex = index;
      }

      tour.add(start);
      Collections.reverse(tour);

      List<String> path = new ArrayList<>();

      for(int i=0; i< tour.size(); i++) {
        path.add(list.get(tour.get(i)).name);
      }

      allPaths.add(path);
      tour = new ArrayList<>();

    }

    System.out.println(allPaths);

  }

  //This method generates all bit sets of size n where r bits 
  // are set to one. The result is returned as a list of integer masks.
  public static List<Integer> combinations(int r, int n) {
    List<Integer> subsets = new ArrayList<>();
    combinations(0, 0, r, n, subsets);
    return subsets;
  }

  //To find all the combinations of size r we need to recurse until we have
  // selected r elements (aka r = 0), otherwise if r != 0 then we still need to select
  // an element which is found after the position of our last selected element
  private static void combinations(int set, int at, int r, int n, List<Integer> subsets) {

    // Return early if there are more elements left to select than what is available.
    int elementsLeftToPick = n - at;
    if (elementsLeftToPick < r) {
        return;
    }

    // We selected 'r' elements so we found a valid subset!
    if (r == 0) {
      subsets.add(set);
    } else {
      for (int i = at; i < n; i++) {
        // Try including this element
        set |= 1 << i;

        combinations(set, i + 1, r - 1, n, subsets);

        // Backtrack and try the instance where we did not include this element
        set &= ~(1 << i);
      }
    }
  }

  private static boolean notIn(int elem, int subset) {
    return ((1 << elem) & subset) == 0;
  }

  private static List<Vertex> reiterate(Map<String, List<Vertex>> group, List<Vertex> vertices) {
    List<Vertex> res = new ArrayList<>();
    Set<String> centroids = group.keySet();

    List<String> list = new ArrayList<String>(centroids);

    double sumx = 0;
    double sumy = 0;

    for(int i=0;i<centroids.size();i++) {
      List<Vertex> current = group.get(list.get(i));
      int n = current.size();
      for(int j=0;j<n;j++){
        sumx+=current.get(j).x;
        sumy+=current.get(j).y;
      }
      res.add(new Vertex(list.get(i), sumx/n, sumy/n));
    }

    return res;
  }

  private static List<Vertex> createCentroids(int clusters, List<Vertex> vertices) {
    List<Vertex> centroids = new ArrayList<>();

    List<Double> xCoords = new ArrayList<>();
    List<Double> yCoords = new ArrayList<>();

    for(Vertex v : vertices) {
      xCoords.add(v.x);
      yCoords.add(v.y);
    }

    Collections.sort(xCoords);
    Collections.sort(yCoords);

    int n = xCoords.size();

    int count = 1;

    int factor = n-1/clusters-1;

    while(count <= clusters) {
      int current = 0;
      centroids.add(new Vertex("v"+count++, xCoords.get(current), yCoords.get(current)));
      current += factor;
    }

    System.out.println(centroids);

    return centroids;
  }

  private static Map<String, List<Vertex>> cluster(List<Vertex> vertices, int k, List<Vertex> centroids, int perCar) {
    PriorityQueue<Edge> pq = new PriorityQueue<Edge>((a,b)-> Double.compare(a.distance, b.distance));
    Set<String> grouped = new HashSet<>();
    Map<String, List<Vertex>> map = new HashMap<>();
    for(Vertex centroid : centroids) {
      map.put(centroid.name, new ArrayList<>());
      for(Vertex vertex : vertices) {
        pq.add(new Edge(centroid, vertex, getDistance(centroid.x, centroid.y, vertex.x, vertex.y)));
      }
    }

    while(!pq.isEmpty()) {
      if(grouped.size() == vertices.size()) {
        break;
      }
      Edge current = pq.poll();
      if(!grouped.contains(current.destination.name) && map.get(current.source.name).size() < perCar) {
        map.get(current.source.name).add(current.destination);
        grouped.add(current.destination.name);
      }
    }

    return map;

  }

  private static double getDistance(double x1, double y1, double x2, double y2) {
    return Point2D.distance(x1, y1, x2, y2);
  }
}