class Vertex
{
    public name:string;
    public x:number;
    public y:number;
    constructor(name: string, x : number, y : number)
    {
        this.name = name;
        this.x = x;
        this.y = y;
    }
    public toString() : string
    {
        return this.name;
    }
    public equals(v : any) : boolean
    {
        if (!(v instanceof Vertex))
        {
            return false;
        }
        return this.x == v.x && this.y==v.y && this.name == v.name;
    }
}

class Edge {
    source:Vertex;
    destination:Vertex;
    distance:number;
  
    constructor(source:Vertex, destination:Vertex, distance:number) {
      this.source = source;
      this.destination = destination;
      this.distance = distance;
    }
  }

class Kmeans
{
    public static  main(args:string[])
    {
        console.log("Hello, World!");
        var v2 = new Vertex("a", 1, 1);
        var v3 = new Vertex("b", 2, 3);
        var v4 = new Vertex("c", -2, 2);
        var v5 = new Vertex("d", -1, 4);
        var v6 = new Vertex("e", -1, -1);
        var v7 = new Vertex("f", 2, -3);
        var v8 = new Vertex("g", -2, -1);
        var v9 = new Vertex("h", 5, 0);
        var v10 = new Vertex("i", 1, 1);
        var v11 = new Vertex("j", 2, 3);
        var v12 = new Vertex("k", -2, 2);
        var v13 = new Vertex("l", -1, 4);
        var v14 = new Vertex("m", 1, 1);
        var v15 = new Vertex("n", 2, 3);
        var v16 = new Vertex("o", -2, 2);
        var v17 = new Vertex("p", -1, 4);
        var vertices[];
        var clusters = 4;
        var perCar = 5;
        var group = cluster(vertices, clusters, createCentroids(clusters, vertices), perCar);
        var group1 = algo.Kmeans.cluster(vertices, clusters, algo.Kmeans.reiterate(group, vertices), perCar);
        while (!java.util.Map.equals(group))
        {
            group = group1;
            group1 = algo.Kmeans.cluster(vertices, clusters, algo.Kmeans.reiterate(group, vertices), perCar);
            console.log(java.util.Map.entrySet());
        }
        var s = new Vertex("s", 0, 0);
        var allPaths = [];
        algo.Kmeans.findAllPaths(allPaths, group1, s);
    }
    private static  findAllPaths(allPaths:[], group1:new Map, startingPoint:Vertex)
    {
        var tour = [];
        var minTourCost = Infinity;
        for ( const  list of group1.values())
        {
            var distance:number[][] = Array(list.size() + 1).fill(0.0).map(()=>new Array(list.size() + 1).fill(0.0));
            list.add(0,startingPoint);
            for (var i = 0; i < list.size() - 1; i++)
            {
                for (var j = i; j < list.size(); j++)
                {
                    distance[i][j] = algo.Kmeans.getDistance(list.get(i).x, list.get(i).y, list.get(j).x, list.get(j).y);
                    distance[j][i] = distance[i][j];
                }
            }
            // start
            var n = distance.length;
            var finalState = parseInt(Math.pow(2,n)) - 1;
            var memo:number[][] = Array(n).fill(0.0).map(()=>new Array(parseInt(Math.pow(2,n))).fill(0.0));
            var start = 0;
            // Add all outgoing edges from the starting node to memo table.
            for (var end = 0; end < n; end++)
            {
                if (end == start)
                {
                    continue;
                }
                memo[end][(1 << start) | (1 << end)] = distance[start][end];
            }
            for (var r = 3; r <= n; r++)
            {
                for ( const  subset of algo.Kmeans.combinations(r, n))
                {
                    if (algo.Kmeans.notIn(start, subset))
                    {
                        continue;
                    }
                    for (var next = 0; next < n; next++)
                    {
                        if (next == start || algo.Kmeans.notIn(next, subset))
                        {
                            continue;
                        }
                        var subsetWithoutNext = subset ^ (1 << next);
                        var minDist = Infinity;
                        for (var end = 0; end < n; end++)
                        {
                            if (end == start || end == next || algo.Kmeans.notIn(end, subset))
                            {
                                continue;
                            }
                            var newDistance = memo[end][subsetWithoutNext] + distance[end][next];
                            if (newDistance < minDist)
                            {
                                minDist = newDistance;
                            }
                        }
                        memo[next][subset] = minDist;
                    }
                }
            }
            // Connect tour back to starting node and minimize cost.
            for (var i = 0; i < n; i++)
            {
                if (i == start)
                {
                    continue;
                }
                var tourCost = memo[i][finalState] + distance[i][start];
                if (tourCost < minTourCost)
                {
                    minTourCost = tourCost;
                }
            }
            var lastIndex = start;
            var state = finalState;
            tour.add(start);
            // Reconstruct TSP path from memo table.
            for (var i = 1; i < n; i++)
            {
                var index = -1;
                for (var j = 0; j < n; j++)
                {
                    if (j == start || algo.Kmeans.notIn(j, state))
                    {
                        continue;
                    }
                    if (index == -1)
                    {
                        index = j;
                    }
                    var prevDist = memo[index][state] + distance[index][lastIndex];
                    var newDist = memo[j][state] + distance[j][lastIndex];
                    if (newDist < prevDist)
                    {
                        index = j;
                    }
                }
                tour.add(index);
                state = state ^ (1 << index);
                lastIndex = index;
            }
            tour.add(start);
            Collections.reverse(tour);
            var path = [];
            for (var i = 0; i < tour.size(); i++)
            {
                path.add(list.get(tour.get(i)).name);
            }
            allPaths.add(path);
            tour = [];
        }
        console.log(allPaths);
    }
    // This method generates all bit sets of size n where r bits
    // are set to one. The result is returned as a list of integer masks.
    public static []<number> combinations(r:number, n:number)
    {
        var subsets = [];
        algo.Kmeans.combinations(0, 0, r, n, subsets);
        return subsets;
    }
    // To find all the combinations of size r we need to recurse until we have
    // selected r elements (aka r = 0), otherwise if r != 0 then we still need to select
    // an element which is found after the position of our last selected element
    private static  combinations(set:number, at:number, r:number, n:number, subsets:[])
    {
        // Return early if there are more elements left to select than what is available.
        var elementsLeftToPick = n - at;
        if (elementsLeftToPick < r)
        {
            return;
        }
        // We selected 'r' elements so we found a valid subset!
        if (r == 0)
        {
            subsets.add(set);
        }
        else 
        {
            for (var i = at; i < n; i++)
            {
                // Try including this element
                set |= 1 << i;
                algo.Kmeans.combinations(set, i + 1, r - 1, n, subsets);
                // Backtrack and try the instance where we did not include this element
                set &= ~(1 << i);
            }
        }
    }
    private static boolean notIn(elem:number, subset:number)
    {
        return ((1 << elem) & subset) == 0;
    }
    private static []<Vertex> reiterate(group:new Map, vertices:[])
    {
        var res = [];
        var centroids = group.keys();
        var list = [];
        var sumx = 0;
        var sumy = 0;
        for (var i = 0; i < centroids.size(); i++)
        {
            var current = group.get(list.get(i));
            var n = current.size();
            for (var j = 0; j < n; j++)
            {
                sumx += current.get(j).x;
                sumy += current.get(j).y;
            }
            res.add(new Vertex(list.get(i), sumx / n, sumy / n));
        }
        return res;
    }
    function [] createCentroids(clusters:number, vertices:[])
    {
        var centroids = [];
        var xCoords = [];
        var yCoords = [];
        for ( const  v of vertices)
        {
            xCoords.add(v.x);
            yCoords.add(v.y);
        }
        Collections.sort(xCoords);
        Collections.sort(yCoords);
        var n = xCoords.size();
        var count = 1;
        var factor = n - parseInt(1 / clusters) - 1;
        while (count <= clusters)
        {
            var current = 0;
            centroids.add(new Vertex("v" + count++, xCoords.get(current), yCoords.get(current)));
            current += factor;
        }
        console.log(centroids);
        return centroids;
    }

    function Map cluster(vertices:[], k:number, centroids:[], perCar:number)
    {
        var pq = java.util.PriorityQueue;
        var grouped = new Set;
        var map = new Map;
        for ( const  centroid of centroids)
        {
            map.set(centroid.name,[]);
            for ( const  vertex of vertices)
            {
                pq.add(new Edge(centroid, vertex, algo.Kmeans.getDistance(centroid.x, centroid.y, vertex.x, vertex.y)));
            }
        }
        while (!pq.isEmpty())
        {
            if (grouped.size() == vertices.size())
            {
                break;
            }
            var current = pq.poll();
            if (!grouped.has(current.destination.name) && map.get(current.source.name).size() < perCar)
            {
                map.get(current.source.name).add(current.destination);
                grouped.add(current.destination.name);
            }
        }
        return map;
    }
    private static number getDistance(x1:number, y1:number, x2:number, y2:number)
    {
        return Point2D.distance(x1,y1,x2,y2);
    }
}
Kmeans.main([]);