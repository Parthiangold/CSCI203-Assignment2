import java.io.BufferedReader;
import java.io.FileReader;
import java.io.InputStreamReader;
import java.io.IOException;


public class Main {

    static final double INF = 1e100;
    static final double EPS = 1e-9;

    // Graph representation (filled after reading)
    static int nVertices;
    static int nEdges;
    static int[] labels;           // labels[i] = vertex label (original int label)
    static double[][] coords;      // coords[i][0]=x, coords[i][1]=y
    static int[] outCount;         // number of outgoing edges per vertex (temporary)
    static int[][] adj;            // adj[u] = array of neighbor indices
    static double[][] adjW;        // adjW[u] = matching weights
    static int startLabel, goalLabel;
    static int startIdx, goalIdx;

    public static void main(String[] args) throws Exception {
        BufferedReader cin = new BufferedReader(new InputStreamReader(System.in));
        System.out.print("Enter graph filename: ");
        String filename = cin.readLine().trim();

        readGraphFromFile(filename);

        // Print basic info
        System.out.println("Number of vertices: " + nVertices);
        System.out.println("Number of edges: " + nEdges);
        System.out.println("Start vertex: " + startLabel);
        System.out.println("Goal vertex: " + goalLabel);
        double euclid = euclideanDistance(startIdx, goalIdx);
        System.out.printf("Euclidean distance between start and goal: %.3f\n", euclid);
        System.out.println();

        // Dijkstra
        Result dij = dijkstra();
        System.out.println("Dijkstra:");
        if (!dij.reachable) {
            System.out.println("Shortest path: No path");
            System.out.println("Length: INF");
        } else {
            System.out.print("Shortest path: ");
            printPathLabels(dij.path);
            System.out.printf("Length: %.3f\n", dij.length);
        }
        System.out.println("Expanded nodes: " + dij.expanded);
        System.out.println();

        // A*
        Result ast = aStar();
        System.out.println("A*:");
        if (!ast.reachable) {
            System.out.println("Shortest path: No path");
            System.out.println("Length: INF");
        } else {
            System.out.print("Shortest path: ");
            printPathLabels(ast.path);
            System.out.printf("Length: %.3f\n", ast.length);
        }
        System.out.println("Expanded nodes: " + ast.expanded);
        System.out.println();

        // DFS longest path
        Result dfsr = dfsLongestPath();
        System.out.println("DFS (Longest Path):");
        if (!dfsr.reachable) {
            System.out.println("Longest path: No path");
            System.out.println("Length: INF");
        } else {
            System.out.print("Longest path: ");
            printPathLabels(dfsr.path);
            System.out.printf("Length: %.3f\n", dfsr.length);
        }
    }

    // ---------------------------
    // Graph input & utilities
    // ---------------------------
    static void readGraphFromFile(String filename) throws IOException {
        BufferedReader br = new BufferedReader(new FileReader(filename));
        String line;

        // Read first non-empty line
        while ((line = br.readLine()) != null && line.trim().isEmpty()) {}
        if (line == null) throw new IOException("Empty file or bad format.");

        String[] parts = line.trim().split("\\s+");
        if (parts.length < 2) throw new IOException("First line must contain nVertices and nEdges.");
        nVertices = Integer.parseInt(parts[0]);
        nEdges = Integer.parseInt(parts[1]);

        labels = new int[nVertices];
        coords = new double[nVertices][2];

        // Read nVertices lines (label x y)
        for (int i = 0; i < nVertices; i++) {
            line = br.readLine();
            while (line != null && line.trim().isEmpty()) line = br.readLine();
            if (line == null) throw new IOException("Unexpected EOF while reading vertices.");
            parts = line.trim().split("\\s+");
            if (parts.length < 3) throw new IOException("Vertex line must have label x y.");
            labels[i] = Integer.parseInt(parts[0]);
            coords[i][0] = Double.parseDouble(parts[1]);
            coords[i][1] = Double.parseDouble(parts[2]);
        }

        // Temporary storage of edges as label triples (we can't assume labels -> indices are consecutive)
        int[] eFromLabel = new int[nEdges];
        int[] eToLabel = new int[nEdges];
        double[] eW = new double[nEdges];

        for (int e = 0; e < nEdges; e++) {
            line = br.readLine();
            while (line != null && line.trim().isEmpty()) line = br.readLine();
            if (line == null) throw new IOException("Unexpected EOF while reading edges.");
            parts = line.trim().split("\\s+");
            if (parts.length < 3) throw new IOException("Edge line must have i j w.");
            eFromLabel[e] = Integer.parseInt(parts[0]);
            eToLabel[e] = Integer.parseInt(parts[1]);
            eW[e] = Double.parseDouble(parts[2]);
            if (eW[e] < 0) throw new IOException("Edge weight must be non-negative.");
        }

        // Last line: start goal
        line = br.readLine();
        while (line != null && line.trim().isEmpty()) line = br.readLine();
        if (line == null) throw new IOException("Missing start/goal line.");
        parts = line.trim().split("\\s+");
        if (parts.length < 2) throw new IOException("Start/Goal line must contain two integers.");
        startLabel = Integer.parseInt(parts[0]);
        goalLabel = Integer.parseInt(parts[1]);

        // Build mapping label -> index via linear search (labels[]).
        startIdx = indexOfLabel(startLabel);
        goalIdx = indexOfLabel(goalLabel);
        if (startIdx == -1 || goalIdx == -1) {
            throw new IOException("Start or goal label not found among vertices.");
        }

        // Count outgoing edges per vertex
        outCount = new int[nVertices];
        for (int e = 0; e < nEdges; e++) {
            int u = indexOfLabel(eFromLabel[e]);
            if (u == -1) {
                throw new IOException("Edge from label " + eFromLabel[e] + " not found.");
            }
            outCount[u]++;
        }

        // Allocate adjacency arrays
        adj = new int[nVertices][];
        adjW = new double[nVertices][];
        for (int i = 0; i < nVertices; i++) {
            adj[i] = new int[outCount[i]];
            adjW[i] = new double[outCount[i]];
            // we'll use a temporary fill pointer
        }

        // fill pointers
        int[] fillPtr = new int[nVertices];
        for (int e = 0; e < nEdges; e++) {
            int u = indexOfLabel(eFromLabel[e]);
            int v = indexOfLabel(eToLabel[e]);
            if (v == -1) {
                throw new IOException("Edge to label " + eToLabel[e] + " not found.");
            }
            int idx = fillPtr[u]++;
            adj[u][idx] = v;
            adjW[u][idx] = eW[e];
        }

        br.close();
    }

    static int indexOfLabel(int lab) {
        for (int i = 0; i < nVertices; i++) {
            if (labels[i] == lab) return i;
        }
        return -1;
    }

    static double euclideanDistance(int uIdx, int vIdx) {
        double dx = coords[uIdx][0] - coords[vIdx][0];
        double dy = coords[uIdx][1] - coords[vIdx][1];
        return Math.sqrt(dx*dx + dy*dy);
    }

    static void printPathLabels(int[] path) {
        if (path == null || path.length == 0) {
            System.out.println("No path");
            return;
        }
        for (int i = 0; i < path.length; i++) {
            System.out.print(path[i]);
            if (i < path.length - 1) System.out.print(" -> ");
        }
        System.out.println();
    }

    // ---------------------------
    // Dijkstra Implementation
    // ---------------------------
    static class Result {
        boolean reachable;
        int[] path;      // sequence of labels
        double length;
        int expanded;
        Result(boolean r, int[] p, double l, int e) {
            reachable = r; path = p; length = l; expanded = e;
        }
    }

    static Result dijkstra() {
        double[] dist = new double[nVertices];
        boolean[] visited = new boolean[nVertices];
        int[] parent = new int[nVertices];
        for (int i = 0; i < nVertices; i++) {
            dist[i] = INF;
            parent[i] = -1;
            visited[i] = false;
        }
        dist[startIdx] = 0.0;
        parent[startIdx] = -1;

        int expanded = 0;

        while (true) {
            int u = minIndexUnvisited(dist, visited);
            if (u == -1) break; // no reachable unvisited node
            if (dist[u] >= INF/2) break;
            // process u
            visited[u] = true;
            expanded++;
            if (u == goalIdx) break; // we can stop when we extract goal in Dijkstra

            // relax neighbors
            for (int k = 0; k < adj[u].length; k++) {
                int v = adj[u][k];
                double w = adjW[u][k];
                double cand = dist[u] + w;
                if (cand + EPS < dist[v]) {
                    dist[v] = cand;
                    parent[v] = u;
                } else if (Math.abs(cand - dist[v]) <= 1e-6) {
                    // tie: pick lexicographically smaller path (by label sequence)
                    // compare candidate path via u vs existing path via parent[v]
                    if (isCandidatePathLexicoSmaller(parent, u, v)) {
                        parent[v] = u;
                    }
                }
            }
        }

        if (dist[goalIdx] >= INF/2) {
            return new Result(false, null, INF, expanded);
        } else {
            int[] seqIdx = reconstructPathIndices(parent, goalIdx);
            int[] seqLabels = indicesToLabels(seqIdx);
            return new Result(true, seqLabels, dist[goalIdx], expanded);
        }
    }

    // choose unvisited vertex with smallest dist
    static int minIndexUnvisited(double[] dist, boolean[] visited) {
        double best = INF;
        int bi = -1;
        for (int i = 0; i < nVertices; i++) {
            if (!visited[i] && dist[i] < best - EPS) {
                best = dist[i];
                bi = i;
            } else if (!visited[i] && Math.abs(dist[i] - best) <= 1e-9 && bi != -1) {
                // tie-breaker by lexicographic on labels? Not strictly necessary here
                if (labels[i] < labels[bi]) bi = i;
            }
        }
        return bi;
    }

    // when candidate path to v (via u) has same length, see if its label sequence is lexicographically smaller
    static boolean isCandidatePathLexicoSmaller(int[] parent, int u, int v) {
        // candidate path indices = pathFromStart... u, v
        int[] cand = reconstructPathIndicesWithExtra(parent, u, v); // returns indices sequence
        int[] curr = reconstructPathIndices(parent, v);
        if (curr == null) return true; // no path currently
        return lexicographicCompareIndices(cand, curr) < 0;
    }

    // compare two index sequences lexicographically by their labels
    // returns -1 if a<b, 0 if equal, +1 if a>b
    static int lexicographicCompareIndices(int[] a, int[] b) {
        int la = a.length, lb = b.length;
        int lm = Math.min(la, lb);
        for (int i = 0; i < lm; i++) {
            int la_lbl = labels[a[i]];
            int lb_lbl = labels[b[i]];
            if (la_lbl < lb_lbl) return -1;
            if (la_lbl > lb_lbl) return 1;
        }
        if (la < lb) return -1;
        if (la > lb) return 1;
        return 0;
    }

    // reconstruct path indices from startIdx to target using parent array (returns indices sequence)
    static int[] reconstructPathIndices(int[] parent, int target) {
        if (parent[target] == -1 && target != startIdx) {
            // if target isn't start and has no parent => unreachable
            return null;
        }
        // collect backwards
        int[] tmp = new int[nVertices];
        int len = 0;
        int cur = target;
        while (cur != -1) {
            tmp[len++] = cur;
            if (cur == startIdx) break;
            cur = parent[cur];
        }
        if (tmp[len - 1] != startIdx) {
            // no path
            return null;
        }
        // reverse
        int[] res = new int[len];
        for (int i = 0; i < len; i++) res[i] = tmp[len - 1 - i];
        return res;
    }

    // reconstruct candidate path indices for v via u (i.e., path to u then v)
    static int[] reconstructPathIndicesWithExtra(int[] parent, int u, int v) {
        if (u == -1) return new int[] { v }; // shouldn't happen
        int[] pathToU = reconstructPathIndices(parent, u);
        if (pathToU == null) {
            return null; // no path to u
        }
        int[] res = new int[pathToU.length + 1];
        for (int i = 0; i < pathToU.length; i++) res[i] = pathToU[i];
        res[res.length - 1] = v;
        return res;
    }

    static int[] indicesToLabels(int[] idxs) {
        if (idxs == null) return null;
        int[] out = new int[idxs.length];
        for (int i = 0; i < idxs.length; i++) out[i] = labels[idxs[i]];
        return out;
    }

    // ---------------------------
    // A* Implementation
    // ---------------------------
    static Result aStar() {
        double[] g = new double[nVertices];
        double[] f = new double[nVertices];
        boolean[] visited = new boolean[nVertices];
        int[] parent = new int[nVertices];

        for (int i = 0; i < nVertices; i++) {
            g[i] = INF;
            f[i] = INF;
            parent[i] = -1;
            visited[i] = false;
        }
        g[startIdx] = 0.0;
        f[startIdx] = g[startIdx] + euclideanDistance(startIdx, goalIdx);
        parent[startIdx] = -1;

        int expanded = 0;

        while (true) {
            int u = minIndexUnvisitedAStar(f, visited, g);
            if (u == -1) break;
            if (g[u] >= INF/2) break;
            visited[u] = true;
            expanded++;
            if (u == goalIdx) break;

            for (int k = 0; k < adj[u].length; k++) {
                int v = adj[u][k];
                double w = adjW[u][k];
                double candG = g[u] + w;
                double candF = candG + euclideanDistance(v, goalIdx);
                if (candG + EPS < g[v]) {
                    g[v] = candG;
                    f[v] = candF;
                    parent[v] = u;
                } else if (Math.abs(candG - g[v]) <= 1e-6) {
                    // equal g: tie-break by lexicographic path sequence to v
                    if (isCandidatePathLexicoSmaller(parent, u, v)) {
                        parent[v] = u;
                        // update f too
                        f[v] = candF;
                    }
                } else if (Math.abs(candF - f[v]) <= 1e-9) {
                    // if f ties, we can consider lexicographic tie-break
                    if (isCandidatePathLexicoSmaller(parent, u, v)) {
                        parent[v] = u;
                        g[v] = candG;
                        f[v] = candF;
                    }
                }
            }
        }

        if (g[goalIdx] >= INF/2) {
            return new Result(false, null, INF, expanded);
        } else {
            int[] seqIdx = reconstructPathIndices(parent, goalIdx);
            int[] seqLabels = indicesToLabels(seqIdx);
            return new Result(true, seqLabels, g[goalIdx], expanded);
        }
    }

    // select unvisited vertex with smallest f; tie-break with smaller g; then label
    static int minIndexUnvisitedAStar(double[] f, boolean[] visited, double[] g) {
        double bestF = INF;
        double bestG = INF;
        int bi = -1;
        for (int i = 0; i < nVertices; i++) {
            if (visited[i]) continue;
            if (f[i] < bestF - EPS) {
                bestF = f[i];
                bestG = g[i];
                bi = i;
            } else if (Math.abs(f[i] - bestF) <= 1e-9) {
                if (g[i] < bestG - EPS) {
                    bestG = g[i];
                    bi = i;
                } else if (Math.abs(g[i] - bestG) <= 1e-9 && bi != -1) {
                    if (labels[i] < labels[bi]) bi = i;
                }
            }
        }
        return bi;
    }

    // ---------------------------
    // DFS for Longest Simple Path
    // ---------------------------
    static Result dfsLongestPath() {
        boolean[] visited = new boolean[nVertices];
        int[] curPathIdx = new int[nVertices]; // indices
        int[] bestPathIdx = new int[nVertices];
        int[] bestLenHolder = new int[1];
        double[] bestWeight = new double[] { -INF }; // store best length
        int[] curLenHolder = new int[1];
        // We'll use recursion
        for (int i = 0; i < nVertices; i++) visited[i] = false;
        curLenHolder[0] = 0;
        bestLenHolder[0] = 0;

        dfsRec(startIdx, visited, curPathIdx, 0, 0.0, bestPathIdx, bestLenHolder, bestWeight);

        if (bestWeight[0] < -INF/2) {
            // no path found
            return new Result(false, null, INF, 0);
        } else {
            int bestLen = bestLenHolder[0];
            int[] seqIdx = new int[bestLen];
            for (int i = 0; i < bestLen; i++) seqIdx[i] = bestPathIdx[i];
            int[] seqLabels = indicesToLabels(seqIdx);
            return new Result(true, seqLabels, bestWeight[0], 0);
        }
    }

    // We'll keep expanded count not required for DFS per spec, but could be added.
    static void dfsRec(int cur, boolean[] visited, int[] curPathIdx, int curLen,
                       double curWeight, int[] bestPathIdx, int[] bestLenHolder, double[] bestWeight) {
        visited[cur] = true;
        curPathIdx[curLen] = cur;
        curLen++;

        if (cur == goalIdx) {
            // reached goal — evaluate
            boolean better = false;
            if (curWeight > bestWeight[0] + 1e-9) {
                better = true;
            } else if (Math.abs(curWeight - bestWeight[0]) <= 1e-9) {
                // tie: lexicographic smaller sequence (labels)
                int[] cand = new int[curLen];
                int[] currBest = new int[bestLenHolder[0]];
                for (int i = 0; i < curLen; i++) cand[i] = curPathIdx[i];
                for (int i = 0; i < bestLenHolder[0]; i++) currBest[i] = bestPathIdx[i];
                if (bestLenHolder[0] == 0) {
                    better = true;
                } else {
                    int cmp = lexicographicCompareIndices(cand, currBest);
                    if (cmp < 0) better = true;
                }
            }
            if (better) {
                bestWeight[0] = curWeight;
                bestLenHolder[0] = curLen;
                // copy curPathIdx into bestPathIdx
                for (int i = 0; i < curLen; i++) bestPathIdx[i] = curPathIdx[i];
            }
            // continue exploring? For longest path we still allow longer paths that revisit goal only as endpoint,
            // but specification says path must be simple and end at goal, so we don't continue from goal.
        } else {
            // continue exploring neighbors
            for (int k = 0; k < adj[cur].length; k++) {
                int v = adj[cur][k];
                double w = adjW[cur][k];
                if (!visited[v]) {
                    dfsRec(v, visited, curPathIdx, curLen, curWeight + w, bestPathIdx, bestLenHolder, bestWeight);
                }
            }
        }

        // backtrack
        visited[cur] = false;
    }
}
