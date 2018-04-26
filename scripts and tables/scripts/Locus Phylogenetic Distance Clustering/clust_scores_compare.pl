#clust_scores_compare.pl - compare the number of clusters generated for different quality scores used:
#1: No clustering, phylogenetic distance cutoff = 0: the merged sequences used to generate phylogenetic trees (CD-hit 90% sequence identity)
#Output the phylogenetic distance cutoff and number of clusters which yields the maximum overall cluster quality score for different calculations:
#2: Q1: Silhouette + Dunn
#3: Q2: Prop. Sequences Clustered + Silhouette + Dunn
#4: Q3: Prop. Sequences Clustered * (Silhouette + Dunn)
#Constrain clustering within the window of Prop. Sequences Clustered 0.1 - 0.9 

use strict;
use warnings;
use Statistics::Basic qw(:all);

my $stat_file;
my @stat_files = <"e:/eps_latest/phylogenetic_trees/clustering/*/*_clusters_stats.txt">;
#my @stat_files = <"./mcl/*/*evodist_clusters_stats.txt">;


my($sys, $sub);

my %dists;
my %max_dists;
my $max_dist;
my $prot1;
my $prot2;

my $dist;
my $boot;


my $dist_cut;
my $clust_num;
my $num_clusts;
my $num_seqs;
my $i;

my %clusts;
my @liner;
my @dists;
my @clusts;
my @prot;
my @d;
my $j;
my $k;

my $prop_clust;
my $sil;
my $dunn;

my %clust_stat;
my $order;

##Step 2: Extract Prop. Sequences Clustered and Avg. Cluster Silhouette Score (Normalized Blast Bitscores) from cluster stats files

foreach $stat_file (@stat_files){
	($sys, $sub) = ($1, $2) if $stat_file =~ /clustering\/(.*)\/(.*)_clusters_stats/;

	$order = 1;
	open STAT, $stat_file;

	while(<STAT>){
		chomp;
		if($. != 1){
			@liner = split /\t/;
			($dist_cut, $num_seqs, $num_clusts, $sil, $dunn) = (@liner[0..2], $liner[8], $liner[11]);
			$prop_clust = ($num_seqs-$num_clusts)/$num_seqs;
			$prop_clust = 1 if $num_clusts == 1;
			
			$clust_stat{$sys}{$sub}{$order} = [$dist_cut, $prop_clust, $sil, $dunn, $num_clusts];
			$order++;
		}
	}
	
	close STAT;
}

##Step 3: Calculate the Overall Clustering Quality Scores


my $outfile = "./cluster_quality_measures.txt";

my $q1;
my $q2;
my $q3;
my %scores;




my @q = ("No_Clust", "Q1: Sil + Dunn", "Q2: Prop_Clust + Sil + Dunn", "Q3: Prop_Clust * (Sil + Dunn)");



#print "System\tSubunit\t",  join("\t", @q), "\n";

#foreach $sys (sort {$a cmp $b} keys %clust_stat){
foreach $sys ("cellulose"){
	foreach $sub (sort {$a cmp $b} keys %{$clust_stat{$sys}}){
		$outfile = "./$sub\_clust_compare.txt";
		
		print $outfile, "\n";
		open OUT, ">", $outfile;
		print OUT join ("\t", @q[1..3]), "\n";
		
		%{$scores{$sys}{$sub}} = ($q[0] => [0, 0, 0],
								$q[1] => [0, 0, 0],
								$q[2] => [0, 0, 0],
								$q[3] => [0, 0, 0]
		);
				
		foreach $order (sort {$a <=> $b} keys %{$clust_stat{$sys}{$sub}}){
			($dist_cut, $prop_clust, $sil, $dunn, $num_clusts) = @{$clust_stat{$sys}{$sub}{$order}};
			
			$q1 = $sil + $dunn;
			$q2 = $sil + $dunn + $prop_clust;
			$q3 = $prop_clust * ($sil + $dunn);
			
			print OUT $q1, "\t", $q2, "\t", $q3, "\n";
			
			if($dist_cut == 0){
				$scores{$sys}{$sub}{$q[0]} = [0, 0, $num_clusts];
			}
			elsif($prop_clust > 0.1 && $prop_clust < 0.9){
				$scores{$sys}{$sub}{$q[1]} = [$dist_cut, $q1, $num_clusts] if $q1 > $scores{$sys}{$sub}{$q[1]}[1];
				$scores{$sys}{$sub}{$q[2]} = [$dist_cut, $q2, $num_clusts] if $q2 > $scores{$sys}{$sub}{$q[2]}[1];
				$scores{$sys}{$sub}{$q[3]} = [$dist_cut, $q3, $num_clusts] if $q3 > $scores{$sys}{$sub}{$q[3]}[1];
			}
		}
		
		close OUT;
		
		#print $sys, "\t", $sub;
		foreach $i (@q){
			#print "\t", $scores{$sys}{$sub}{$i}[0];
		}
		#print "\n";
	}
}


