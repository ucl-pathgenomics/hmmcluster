
package hmmercluster;

// New August model from Richard. model1 = BIC, model2 - bic, TestEntropyModel = this
// This model was not used in the PNAS publication, but kept here for interest or future endeavours

import java.io.FileReader;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.HashMap;
import java.util.TreeMap;
import java.util.TreeSet;
import pal.alignment.Alignment;
import pal.alignment.AlignmentReaders;
import pal.datatype.Nucleotides;

public class TestEntropyModel {

    Alignment align = null;                                         // Alignment
    String[] identifiers = null;                                    // Identifiers
    int[][] sequences = null;                                       // Sequences in int format
    int nSeqs = 0;                                                  // Number of sequences
    int nSites = 0;                                                 // Number of sites
    HashMap<Character, Integer> baseToInt = new HashMap<>();        // Mapping of bases to ints, gap = 4, other = 5
    TreeMap<String, Cluster> computedClusterMap = new TreeMap<>();  // Repository for previously computed values
    boolean verbose = false;                                        // Output
    boolean trim = false;                                            // Trim ends


    double maxFracPoly = 0.6;
    double minScaledEntropy = -0.3;
    double minOlap = 0.6;
    int iStartLimit = 4;
    int iEndLimit = 1024;

    public static void main(String[] args) {
        TestEntropyModel hmm = new TestEntropyModel(args);
        hmm.run();
    }

    TestEntropyModel(String[] args) {
        readFiles(args);                                            // Read in alignment
    }

    void run() {
        double best = 10.0;
        int iStart = 0;
        int iEnd = nSites - 1;
        int iStartBest = 0;
        int iEndBest =  nSites;

        System.out.println(iEnd);
        if (!trim) {
            iStartBest = iStartLimit;
            iEndBest = nSites;
        } else {
            int iStartBestCoarse = iStartBest;
            int iEndBestCoarse = iEndBest;
            double bestCoarse = 10.0;
            for (int iStartTry = 0; iStartTry < iEnd-5; iStartTry+=5) {
                double newValue = findBest(iStartTry, iEnd, false);
                System.out.println(iStartTry + "\t" + iEnd + "\t" + newValue);
                if (newValue < bestCoarse) {
                    bestCoarse = newValue;
                    iStartBestCoarse = iStartTry;
                }
            }
            for (int iEndTry = iStartBest+5; iEndTry < iEnd; iEndTry+=5) {
                double newValue = findBest(iStartBest, iEndTry, false);
                System.out.println(iStartBest + "\t" + iEndTry + "\t" + newValue);
                if (newValue < bestCoarse) {
                    bestCoarse = newValue;
                    iEndBestCoarse = iEndTry;
                }
            }
            for (int iStartTry = Math.max(0, iStartBestCoarse-20); iStartTry < Math.min(iEndBestCoarse-5, iStartBestCoarse+20); iStartTry+=1) {
                double newValue = findBest(iStartTry, iEnd, false);
                System.out.println(iStartTry + "\t" + iEnd + "\t" + newValue);
                if (newValue < best) {
                    best = newValue;
                    iStartBest = iStartTry;
                }
            }
            for (int iEndTry = Math.max(iStartBest+5, iEndBestCoarse-20); iEndTry < Math.min(iEnd, iEndBestCoarse+20); iEndTry+=1) {
                double newValue = findBest(iStartBest, iEndTry, false);
                System.out.println(iStartBest + "\t" + iEndTry + "\t" + newValue);
                if (newValue < best) {
                    best = newValue;
                    iEndBest = iEndTry;
                }
            }
        }

        findBest(iStartBest, iEndBest, true);
    }

    /**
     *
     * @param iStart Start of region
     * @param iEnd End of region
     * @param printOutput Print out stuff
     * @return  Best BIC value relative to single cluster
     */
    double findBest(int iStart, int iEnd, boolean printOutput) {

        TreeSet<Integer> all = new TreeSet<>();                     // Single cluster
        ArrayList<Cluster>[] current = new ArrayList[nSeqs+1];      // Repository for sets of clusters
        computedClusterMap.clear();                                 // Scrub memory

        current[0] = new ArrayList<>();                             // Initial set of clusters
        for (int iSeq = 0; iSeq < nSeqs; iSeq++) {
            TreeSet<Integer> example = new TreeSet<>();
            example.add(iSeq);
            String label = Arrays.toString(example.toArray()).replaceAll(" ", "");
            Cluster newCluster = new Cluster(example, label, iStart, iEnd);
            current[0].add(newCluster);                             // Add individual sequences to current set
            computedClusterMap.put(label, newCluster);              // Add to memory
            all.add(iSeq);                                          // Add individual sequences to all cluster
        }
        String label = Arrays.toString(all.toArray()).replaceAll(" ", "");
        Cluster allCluster = new Cluster(all, label, iStart, iEnd);
        computedClusterMap.put(label, allCluster);

        // Which iteration yields best results
        int iIter;
        for (iIter = 1; iIter < sequences.length; iIter++) {                            // Iterate combining sets
            current[iIter] = iterate(current[iIter-1], iStart, iEnd);
            if (current[iIter] == null) {
                break;
            }

        }


        double totalScaledEntropy = 0.0;
        for (Cluster cluster : current[iIter-1]) {           // Sum quantities over clusters
            totalScaledEntropy += cluster.scaledEntropy;
        }
        if (verbose) {
            System.out.println("ggg\t" + current[iIter-1].size() + "\t" +
                    totalScaledEntropy + "\t" + (totalScaledEntropy/current[iIter-1].size()) + "\t" + allCluster.scaledEntropy);
        }

        if (printOutput) {
            System.out.format("%d\t%d\t%d\t%.6f\t%.6f\n", current[iIter-1].size(), iStart, iEnd,
                    totalScaledEntropy, (totalScaledEntropy - allCluster.scaledEntropy));
            for (Cluster cluster : current[iIter-1]) {
                System.out.format("\n%s\t%.6f\t%.6f\t%.0f\t%.2f\n", cluster.label, cluster.entropy, cluster.scaledEntropy,
                        cluster.nActiveSites, cluster.fracPoly);
                for (int iSeq : cluster.includedSequences) {
                    System.out.println(iSeq + "\t" + identifiers[iSeq]);
                }
            }
        }


        return totalScaledEntropy/current[iIter-1].size()-allCluster.scaledEntropy;
    }



    /**
     * Compute change in BIC compared with single cluster and print stuff
     * @param current               // Set of clusters
     * @param allCluster            // Single cluster
     * @param printOutput           // Print stuff
     * @param printClusters         // Print identities of clusters
     * @return change in BIC compared with single cluster
     */
//    double summariseCurrent(ArrayList<Cluster> current, Cluster allCluster, int iStart, int iEnd, boolean printOutput, boolean printClusters) {
//        double totalEntropy = 0.0;
//        double allEntropy = allCluster.entropy;
//        for (Cluster cluster : current) {           // Sum quantities over clusters
//            totalEntropy += cluster.entropy;
//        }
//        
//        if (printOutput) {
//            if (printClusters) {
//                System.out.println();
//            }
//            System.out.format("%d\t%d\t%d\t%.2f\t%.0f\t%.2f\t%.2f\n", current.size(), iStart, iEnd, totalEntropy, totalParams, totalBIC, (totalBIC - allEntropy));
//            if (printClusters) {
//                for (Cluster cluster : current) {
//                    System.out.format("\n%s\t%.2f\t%.0f\t%.2f\n", cluster.label, cluster.entropy, cluster.nParams, cluster.score);
//                    for (int iSeq : cluster.includedSequences) {
//                        System.out.println(iSeq + "\t" + identifiers[iSeq]);
//                    }
//                }
//            }
//        }
//        return (totalBIC - allEntropy);
//    }


    /**
     * Find best way to fuse two clusters 
     * @param current Current set of clusters
     * @param iStart Start of aligned region
     * @param iEnd End of aligned region
     * @return New combined cluster
     */
    ArrayList<Cluster> iterate(ArrayList<Cluster> current, int iStart, int iEnd) {
        double best = -9999999999.9999;
        int bestI = -9;
        int bestJ = -9;

        TreeSet<Integer> included = new TreeSet<>();                                // Indices of combined set
        String[][] labels = new String[current.size()][current.size()];             // Labels for all possible combinations of two sets

        for (int iSeq = 0; iSeq < current.size()-1; iSeq++) {                       // Sum over all pairs of sets
            for (int jSeq = iSeq+1; jSeq < current.size(); jSeq++) {
                double previousEntropy = current.get(iSeq).entropy + current.get(jSeq).entropy; // Previous BIC sum for components of combined set
                included.clear();                                                   // Assemble indices of combined sets
                included.addAll(current.get(iSeq).includedSequences);
                included.addAll(current.get(jSeq).includedSequences);
                labels[iSeq][jSeq] = Arrays.toString(included.toArray()).replaceAll(" ", "");
                Cluster newCluster;                                                 // New cluster
                if (computedClusterMap.containsKey(labels[iSeq][jSeq])) {
                    newCluster = computedClusterMap.get(labels[iSeq][jSeq]);
                } else {
                    newCluster = new Cluster(included, labels[iSeq][jSeq], iStart, iEnd);
                    computedClusterMap.put(labels[iSeq][jSeq], newCluster);
                }
                double deltaEntropy = newCluster.entropy - previousEntropy;                     // Change in BIC
                double oLap = ((current.get(iSeq).nActiveSites + current.get(jSeq).nActiveSites
                        - newCluster.nActiveSites)/newCluster.nActiveSites);
//                System.out.println(labels[iSeq][jSeq] + "\t" + newCluster.entropy + "\t" 
//                        + current.get(iSeq).entropy + "\t" + current.get(jSeq).entropy
//                        + "\t" + deltaEntropy + "\t" + newCluster.fracPoly 
//                        + "\t" + newCluster.scaledEntropy
//                        + "\t" + newCluster.nActiveSites + "\t" 
//                        + "\t" + current.get(iSeq).nActiveSites + "\t" + current.get(jSeq).nActiveSites
//                        + "\t" + (current.get(iSeq).nActiveSites + current.get(jSeq).nActiveSites)
//                        + "\t" + ((current.get(iSeq).nActiveSites + current.get(jSeq).nActiveSites
//                                - newCluster.nActiveSites)/newCluster.nActiveSites));
                if (oLap >= minOlap && newCluster.fracPoly <= maxFracPoly &&
                        newCluster.scaledEntropy >= minScaledEntropy &&
                        deltaEntropy > best) {                                              // Best reduction in BIC
                    best = deltaEntropy;
                    bestI = iSeq;
                    bestJ = jSeq;
                }
            }
        }

        if (bestI >= 0 && bestJ >= 0) {
            ArrayList<Cluster> newCurrent = new ArrayList<>();                          // Create new set of clusters
            for (int iSeq = 0; iSeq < current.size(); iSeq++) {                         // Copy over unaffected clusters
                if (iSeq != bestI && iSeq != bestJ) {
                    newCurrent.add(current.get(iSeq));
                }
            }
            Cluster newCluster = computedClusterMap.get(labels[bestI][bestJ]);          // Add new cluster
            newCurrent.add(newCluster);
            if (verbose) {
                System.out.println("kkk\t" + newCurrent.size() + "\t" + newCluster.label + "\t" + newCluster.entropy + "\t" + newCluster.fracPoly
                        + "\t" + newCluster.scaledEntropy);
            }
            return newCurrent;
        }
        return null;                                                // Return new set of clusters
    }

    /**
     * Read in alignment
     * @param args Name of alignment file
     */
    void readFiles(String[] args) {
        // Load up baseToInt
        baseToInt.put('A', 0);
        baseToInt.put('C', 1);
        baseToInt.put('G', 2);
        baseToInt.put('T', 3);
        baseToInt.put('a', 0);
        baseToInt.put('c', 1);
        baseToInt.put('g', 2);
        baseToInt.put('t', 3);
        baseToInt.put('-', 4);

        // Read in alignment
        try {
            Nucleotides dataType = new Nucleotides();
            FileReader in = new FileReader(args[0]);
            align = AlignmentReaders.readFastaSequences(in, dataType);
            in.close();
        } catch (IOException e) {
            System.out.println("Error: File not found (IO error)");
            System.exit(1);
        }
        nSites = align.getSiteCount();
        nSeqs = align.getSequenceCount();
        System.out.println("Number of sites:\t" + nSites);
        System.out.println("Number of sequences:\t" + nSeqs);

        // Fill sequences and identifiers
        sequences = new int[nSeqs][nSites];
        identifiers = new String[nSeqs];
        for (int iSeq = 0; iSeq < nSeqs; iSeq++) {
            identifiers[iSeq] = align.getIdentifier(iSeq).toString();
            for (int iSite = 0; iSite < nSites; iSite++) {
                sequences[iSeq][iSite] = 5;
                if (baseToInt.containsKey(align.getData(iSeq, iSite))) {
                    sequences[iSeq][iSite] = baseToInt.get(align.getData(iSeq, iSite));
                }
            }
        }
    }

    /**
     * Class containing group of similar sequences
     */
    class Cluster {
        TreeSet<Integer> includedSequences = new TreeSet<>();       // Indices of included sequences
        double entropy = 0.0;                                            // log likelihood of HMM representation
        double nActiveSites = 0.0;                                  // Number of parameters                                          // BIC
        String label = "";                                          // String version of includedSequences
        double scaledEntropy = 0.0;
        double fracPoly = 0.0;

        /**
         * Initialise new cluster
         * @param addSequences Indices of sequences to be added
         * @param label Indices in String format
         * @param iStart Start of aligned regions
         * @param iEnd End of aligned regions
         */
        Cluster(TreeSet<Integer> addSequences, String label, int iStart, int iEnd) {
            includedSequences.clear();
            includedSequences.addAll(addSequences);
            this.label = label;
            computeLL(iStart, iEnd);
        }

        void computeLL(int iStart, int iEnd) {
            for (int iSite = iStart; iSite < iEnd; iSite++) {
                int[] baseCount = new int[4];                               // Count of each type of base
                int totBaseCount = 0;                                       // Total count

                for (int iSeq : includedSequences) {
                    if (sequences[iSeq][iSite] < 4) {                       // If standard base
                        baseCount[sequences[iSeq][iSite]]++;
                        totBaseCount++;
                    }
                }

                int nBasePresent = 0;                                       // How many bases are present
                for (int iBase = 0; iBase < 4; iBase++) {
                    if (baseCount[iBase] > 0) {                             // If base is present
                        double freq = baseCount[iBase]*1.0/totBaseCount;    // Freq of base
                        entropy += freq * Math.log(freq);            // Contribution to log likelihood
                        nBasePresent++;
                    }
                }
                if (nBasePresent > 0) {
                    nActiveSites++;
                    if (nBasePresent > 1) {
                        fracPoly++;
                    }
                }
            }
            if (nActiveSites > 0.01) {
                scaledEntropy = entropy / nActiveSites;
                fracPoly /= nActiveSites;
            }
        }
    }

}