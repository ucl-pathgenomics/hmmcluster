
package hmmercluster;

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

/*
Author: Richard Goldstein
Edited: Oscar Charles
July 2020

Clusters regions on CMV genomes using AIC and hmm data.

*/



/* Simple Process Breakdown
class construct
HMMCluster
    readFiles()
        takes alignment file location
        reads the fasta into a character matrix
        returns "sequences" a matrix where atcg- -> 12345
run




********* oscar notes
loops over start of 5 -> end - 5.  given end =  end. gets LL of clusters? for this alignment region?
then
loops over end 0 -> end - 5 given start  = 0. Gets LL of clusters? for this alignment region?


*/

//todo make the searching mechanism faster rather than by +5. do by 50 for 10* speed gain. then more granular?

public class HMMCluster{

    Alignment align = null;                                         // Alignment
    String[] identifiers = null;                                    // Identifiers
    int[][] sequences = null;                                       // Sequences in int format
    int nSeqs = 0;                                                  // Number of sequences
    int nSites = 0;                                                 // Number of sites
    HashMap<Character, Integer> baseToInt = new HashMap<>();        // Mapping of bases to ints, gap = 4, other = 5
    TreeMap<String, Cluster> computedClusterMap = new TreeMap<>();  // Repository for previously computed values
    boolean verbose = false;                                        // Output
    boolean trim = true;                                            // Trim ends

    public static void main(String[] args) {
        HMMCluster hmm = new HMMCluster(args);                      // create a new object of class HMMCluster, then run that object
        hmm.run();
    }
    
    HMMCluster(String[] args) {
        readFiles(args);                                            // Read in alignment
    }

    void run() {
        double best = 10.0;
        int iStart = 0;
        int iEnd = nSites;
        int iStartBest = 0;
        int iEndBest = nSites;
        int iGap = 5;
        System.out.print("Coarse gap:" + "\t" + iGap + "t");
        if (trim) { // ===TRUE
            int iStartBestCoarse = iStartBest;
            int iEndBestCoarse = iEndBest;
            double bestCoarse = 10.0;
            for (int iStartTry = 0; iStartTry < iEnd-iGap; iStartTry+=iGap) {     // from 5 to end - 5 : end
                double newValue = findBest(iStartTry, iEnd, false);
                System.out.println(iStartTry + "\t" + iEnd + "\t" + newValue);
                if (newValue < bestCoarse) {
                    bestCoarse = newValue;
                    iStartBestCoarse = iStartTry;
                }
            }
            for (int iEndTry = iStartBest+iGap; iEndTry < iEnd; iEndTry+=iGap) {  // start : from 5 to end
                double newValue = findBest(iStartBest, iEndTry, false);
                System.out.println(iStartBest + "\t" + iEndTry + "\t" + newValue);
                if (newValue < bestCoarse) {
                    bestCoarse = newValue;
                    iEndBestCoarse = iEndTry;
                }
            }
            for (int iStartTry = Math.max(0, iStartBestCoarse-20); iStartTry < Math.min(iEndBestCoarse-iGap, iStartBestCoarse+20); iStartTry+=1) {
                double newValue = findBest(iStartTry, iEnd, false);
                System.out.println(iStartTry + "\t" + iEnd + "\t" + newValue);
                if (newValue < best) {
                    best = newValue;
                    iStartBest = iStartTry;
                }
            }
            for (int iEndTry = Math.max(iStartBest+iGap, iEndBestCoarse-20); iEndTry < Math.min(iEnd, iEndBestCoarse+20); iEndTry+=1) {
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
     * @return  Best AIC value relative to single cluster
     */
    double findBest(int iStart, int iEnd, boolean printOutput) {
        
        TreeSet<Integer> all = new TreeSet<>();                     // Single cluster <tree class from pal>
        ArrayList<Cluster>[] current = new ArrayList[nSeqs+1];      // Repository for sets of clusters
        computedClusterMap.clear();                                 // Scrub memory
        
        current[0] = new ArrayList<>();                             // Initial set of clusters
        for (int iSeq = 0; iSeq < nSeqs; iSeq++) {
            TreeSet<Integer> example = new TreeSet<>();
            example.add(iSeq);
            String label = Arrays.toString(example.toArray());
            Cluster newCluster = new Cluster(example, label, iStart, iEnd); 
            current[0].add(newCluster);                             // Add individual sequences to current set
            computedClusterMap.put(label, newCluster);              // Add to memory
            all.add(iSeq);                                          // Add individual sequences to all cluster
        }
        String label = Arrays.toString(all.toArray());
        Cluster allCluster = new Cluster(all, label, iStart, iEnd);
        computedClusterMap.put(label, allCluster);
        
        // Find best set of clusters using greedy algorithm
        double bestScore = summariseCurrent(current[0], allCluster, 
                iStart, iEnd, printOutput, verbose);                                        // Start with independent sequences
        int bestIter = 0;                                                                   // Which iteration yields best results
        for (int iIter = 1; iIter < sequences.length; iIter++) {                            // Iterate combining sets
            current[iIter] = iterate(current[iIter-1], iStart, iEnd);
            double newScore = summariseCurrent(current[iIter], allCluster, 
                    iStart, iEnd, printOutput, verbose);                                    // AIC relative to single cluster
            if (newScore < bestScore) {
                bestScore = newScore;
                bestIter = iIter;
            }
        }
        summariseCurrent(current[bestIter], allCluster, 
                iStart, iEnd, printOutput, true);                                           // Summarise current printing results if printOutput
        return bestScore;
    }   
    
    

    /**
     * Compute change in AIC compared with single cluster and print stuff
     * @param current               // Set of clusters
     * @param allCluster            // Single cluster
     * @param printOutput           // Print stuff
     * @param printClusters         // Print identities of clusters
     * @return change in AIC compared with single cluster
     */
    double summariseCurrent(ArrayList<Cluster> current, Cluster allCluster, int iStart, int iEnd, boolean printOutput, boolean printClusters) {
        double totalLL = 0.0;
        double totalParams = 0.0;
        double totalAIC = 0.0;
        double allAIC = allCluster.AIC;
        for (Cluster cluster : current) {           // Sum quantities over clusters
            totalLL += cluster.ll;
            totalParams += cluster.nParams;
            totalAIC += cluster.AIC;
        }
        
        if (printOutput) {
            if (printClusters) {
                System.out.println();
            }
            System.out.format("%d\t%d\t%d\t%.2f\t%.0f\t%.2f\t%.2f\n", current.size(), iStart, iEnd, totalLL, totalParams, totalAIC, (totalAIC - allAIC));
            if (printClusters) {
                for (Cluster cluster : current) {
                    System.out.format("\n%s\t%.2f\t%.0f\t%.2f\n", cluster.label, cluster.ll, cluster.nParams, cluster.AIC);
                    for (int iSeq : cluster.includedSequences) {
                        System.out.println(iSeq + "\t" + identifiers[iSeq]);
                    }
                }
            }
        }
        return (totalAIC - allAIC);
    }
    
    /** 
     * Find best way to fuse two clusters 
     * @param current Current set of clusters
     * @param iStart Start of aligned region
     * @param iEnd End of aligned region
     * @return New combined cluster
     */
    ArrayList<Cluster> iterate(ArrayList<Cluster> current, int iStart, int iEnd) { // figures out all the ways of pushing things together
        double best = 9999999999.9999;
        int bestI = -9;
        int bestJ = -9;

        TreeSet<Integer> included = new TreeSet<>();                                // Indices of combined set
        String[][] labels = new String[current.size()][current.size()];             // Labels for all possible combinations of two sets

        for (int iSeq = 0; iSeq < current.size()-1; iSeq++) {                       // Sum over all pairs of sets
            for (int jSeq = iSeq+1; jSeq < current.size(); jSeq++) {
                double previousAIC = current.get(iSeq).AIC + current.get(jSeq).AIC; // Previous AIC sum for components of combined set
                included.clear();                                                   // Assemble indices of combined sets
                included.addAll(current.get(iSeq).includedSequences);
                included.addAll(current.get(jSeq).includedSequences);
                labels[iSeq][jSeq] = Arrays.toString(included.toArray());
                Cluster newCluster;                                                 // New cluster
                if (computedClusterMap.containsKey(labels[iSeq][jSeq])) {
                    newCluster = computedClusterMap.get(labels[iSeq][jSeq]);
                } else {
                    newCluster = new Cluster(included, labels[iSeq][jSeq], iStart, iEnd);
                    computedClusterMap.put(labels[iSeq][jSeq], newCluster);
                }
                double deltaAIC = newCluster.AIC - previousAIC;                     // Change in AIC
//                System.out.println(labels[iSeq][jSeq] + "\t" + newCluster.ll + "\t" 
//                        + newCluster.nParams + "\t" + newCluster.AIC + "\t" + deltaAIC);

                if (deltaAIC < best) {                                              // Best reduction in AIC
                    best = deltaAIC;
                    bestI = iSeq;
                    bestJ = jSeq;
                }
            }
        }
        
        ArrayList<Cluster> newCurrent = new ArrayList<>();                          // Create new set of clusters
        for (int iSeq = 0; iSeq < current.size(); iSeq++) {                         // Copy over unaffected clusters
            if (iSeq != bestI && iSeq != bestJ) {
                newCurrent.add(current.get(iSeq));
            }
        }
        Cluster newCluster = computedClusterMap.get(labels[bestI][bestJ]);          // Add new cluster
        newCurrent.add(newCluster);
        return newCurrent;                                                          // Return new set of clusters
    }
    
    /**
     * Read in alignment
     * @param args Name of alignment file
     * @return matrix of atcg- -> 12345
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
    class Cluster { // array of integers
        // Treeset if a vector but sorted by value sorted vector. small -> large.
        TreeSet<Integer> includedSequences = new TreeSet<>();       // Indices of included sequences
        double ll = 0.0;                                            // log likelihood of HMM representation
        double nParams = 0.0;                                       // Number of parameters
        double AIC = 0.0;                                           // AIC
        String label = "";                                          // String version of includedSequences

        // given set of sequences,

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
        
        void computeLL(int iStart, int iEnd) { // look at writeup
            for (int iSite = iStart; iSite < iEnd; iSite++) {
                int[] baseCount = new int[4];                               // Count of each type of base
                int totBaseCount = 0;                                       // Total count
                int[] newStartDelCount = new int[2];                        // Match to match, match to delete
                int[] newEndDelCount = new int[2];                          // Delete to delete, delete to match

                /*
                for each site in alignment how many acgt's
                then add to the log likelihood the number of a's * log(freq(a)* + c * ... * t

                 */


                for (int iSeq : includedSequences) {
                    if (sequences[iSeq][iSite] < 4) {                       // If standard base
                        baseCount[sequences[iSeq][iSite]]++;
                        totBaseCount++;
                    }
                    if (iSite == 0 || sequences[iSeq][iSite-1] != 4) {      // Count match to delete, match to match
                        if (sequences[iSeq][iSite] == 4){
                            newStartDelCount[1]++;
                        } else {
                            newStartDelCount[0]++;
                        }
                    }  
                    if (iSite != nSites-1 && sequences[iSeq][iSite] == 4) { // Count delete to match, delete to delete
                        if (sequences[iSeq][iSite+1] != 4){
                            newEndDelCount[1]++;
                        } else {
                            newEndDelCount[0]++;
                        }
                    }     
                }
                
                int nBasePresent = 0;                                       // How many bases are present

                // hi oscar
                for (int iBase = 0; iBase < 4; iBase++) {                
                    if (baseCount[iBase] > 0) {                             // If base is present
                        double freq = baseCount[iBase]*1.0/totBaseCount;    // Freq of base
                        ll += baseCount[iBase] * Math.log(freq);            // Contribution to log likelihood
                        nBasePresent++;
                    }
                    // hi oscar end
                }
                nParams += Math.min(nBasePresent, 3) + 1;                   // One param for each base, plus one for start of delete
                if (newEndDelCount[0] + newEndDelCount[1] > 0) {            // If delete state is relevant, add parameter for return to match
                    nParams += 1.0;
                }
                if (newStartDelCount[0] > 0 && newStartDelCount[1] > 0) {   // Log likelihood for start of gap
                    double frac0 = newStartDelCount[0]*1.0/(newStartDelCount[0]+newStartDelCount[1]);
                    ll += newStartDelCount[0] * Math.log(frac0) + newStartDelCount[1] * Math.log(1.0-frac0);
                }
                 if (newEndDelCount[0] > 0 && newEndDelCount[1] > 0) {      // Log likelihood for end of gap
                    double frac0 = newEndDelCount[0]*1.0/(newEndDelCount[0]+newEndDelCount[1]);
                    ll += newEndDelCount[0] * Math.log(frac0) + newEndDelCount[1] * Math.log(1.0-frac0);
                }
            } 
            nParams ++;                                                     // Extra parameter for probability of entering model
            ll += includedSequences.size() *                                // Probability of sequence model
                    (Math.log(includedSequences.size()) - Math.log(nSeqs));
            AIC = 2.0 * nParams - 2.0 * ll;                                 // AIC
        }  
    }
}