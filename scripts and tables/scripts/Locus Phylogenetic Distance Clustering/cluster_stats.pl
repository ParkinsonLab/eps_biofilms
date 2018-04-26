#cluster_stats.pl - calculate internal cluster quality validation measures (silhouette coefficient, dunn index) and others ()
#for each phylogenetic clustering generated over a range of rate a.a. change/site cutoffs.
#Include weighted sum of squares of mean differences within and between clusters?

use strict;
use warnings;

#my @infiles = <"./*/*_dists.txt">;
my @infiles = <"./clustering_results/*/*_clusters.txt">;

my @op_files = <"../Operon\ Reconstruction/*_op_final.txt">;

my $outdir = "tmp";

my $outfile;

my $clust_file;
my $dist_file;

my($gi, $gi1, $gi2, $dist);
my %dists;

my @liner;

my @dists;
my $dist_cut;
my $clust_num;

my $i;
my $j;
my @clusts;
my %clusts;

my(@gi1, @gi2);

my $num_clusts;
my $num_seqs;

my $intra_dist;
my $inter_dist;
my $avg_intra_dist;
my $avg_inter_dist;
my $dunn;
my $avg_sil_width;
my $simp;
my $recip_simp;
my($intra_shannon1, $inter_shannon1, $max_shannon1);
my($intra_shannon2, $inter_shannon2, $max_shannon2);

my $op_file;
my $tax;
my $tax_cut = 1;
my @tax;
my %tax_dist;

foreach $op_file (@op_files){
	open OP, $op_file;
	
	while(<OP>){
		chomp;
		@liner = split /\t/;
		$tax = $liner[1];
		@gi1 = @liner[4..$#liner];
		
		foreach $gi (@gi1){
			$gi = $1 if $gi =~ /\|(\d*)\|/;
			$tax_dist{$gi} = $tax;
			
			@tax = split /; /, $tax_dist{$gi};
		}
	}
	
	close OP;

}

my $bitsc;
my $tmp_file;

foreach $clust_file (@infiles){
	#if($dist_file !~ /Bcs/){
	$dist_file = $clust_file;
	$dist_file =~ s/_clusters/_dists/;
	$dist_file =~ s/clustering_results/phylogenetic_distances/;
	
	
	$outfile = $clust_file;
	$outfile =~ s/clustering_results.*\//$outdir\//;
	
	$outfile =~ s/_clusters/_clusters_stats/;
	
	
	%dists = ();
	%clusts = ();

	open DISTS, $dist_file or die;
	
	while(<DISTS>){
		chomp;
		if($. != 1){
			($gi1, $gi2, $dist) = split /\t/;
			$dists{$gi1}{$gi2} = $dist;
			$dists{$gi2}{$gi1} = $dist;

		}
	}
	
	close DISTS;
	
	#@gi1 = keys %dists;
	#@gi2 = @gi1;
	#calc_avg_dist(\@gi1, \@gi2, $dist_file);

	$num_seqs = 0;
	$num_clusts = 0;
	
	if(1){
	#print $clust_file, "\n";
	open CLUSTS, $clust_file;
	
	while(<CLUSTS>){
		chomp;
		@liner = split /\t/;
		
		if($. == 1){
			@dists = splice @liner, 1, $#liner;
			#print $dists[0], "\t", $dists[$#dists], "\n";
		}
		else{
			$gi = $liner[0];
			@clusts = @liner[1..$#liner];
			if(exists $tax_dist{$gi}){
				for($i = 0; $i <= $#clusts; $i++){
					$clust_num = $clusts[$i];
					$dist_cut = $dists[$i];
				
					$clusts{$dist_cut}{$clust_num}{$gi} = 1;
				}
			
				$num_seqs++;
			}
		}
	}
	
	close CLUSTS;
	
	#$dist_cut = 0;
	
	#calculate dunn, silhouette and simpson index for a given distance cutoff threshold clustering from minimum inter-cluster and maximum intra-cluster distances
	
	print $outfile, "\n";
	
	open OUT, ">", $outfile;
	
	print OUT "#dist_cut\tnum_seqs\tnum_clusters\tnum_clusters(recip_simpson_index)\tavg_intra_clust_phylo_dist_phylum(norm_shannon_index)\tavg_inter_clust_phylo_dist_phylum(norm_shannon_index)\tavg_intra_clust_phylo_dist_order(norm_shannon_index)\tavg_inter_clust_phylo_dist_order(norm_shannon_index)\tavg_silhouette_width\tavg_inter_clust_dist\tavg_intra_clust_dist\tdunn_index\n";
	
	@clusts = keys %{$clusts{$dists[$#dists]}};
	
	my $tmp;
	($max_shannon1, $tmp) = calc_shannon($dists[$#dists], $#clusts+1, \@clusts, 1);
	($max_shannon2, $tmp) = calc_shannon($dists[$#dists], $#clusts+1, \@clusts, 3);
	$max_shannon1 = 1 if $max_shannon1 == 0;
	$max_shannon2 = 1 if $max_shannon2 == 0;

	foreach $dist_cut (@dists){

		@clusts = sort {$a <=> $b} keys %{$clusts{$dist_cut}};
		$num_clusts = $#clusts + 1;

		($avg_inter_dist, $avg_intra_dist, $dunn) = calc_dunn($dist_cut, \@clusts);
		$avg_sil_width = calc_silhouette($dist_cut, \@clusts);
		$recip_simp = calc_simpson($dist_cut, $num_seqs, $num_clusts);
		($intra_shannon1, $inter_shannon1) = calc_shannon($dist_cut, $num_clusts, \@clusts, 1);
		#if($dist_cut == 0.900000000000001 && $outfile =~ /BcsA/){
		($intra_shannon2, $inter_shannon2) = calc_shannon($dist_cut, $num_clusts, \@clusts, 3);
		#}
		#print $max_shannon, "\n";
		print OUT $dist_cut, "\t", $num_seqs, "\t", $num_clusts, "\t", $recip_simp, "\t", 
			$intra_shannon1/$max_shannon1, "\t", $inter_shannon1/$max_shannon1, "\t", $intra_shannon2/$max_shannon2, "\t", $inter_shannon2/$max_shannon2, "\t",
			$avg_sil_width, "\t", 
			$avg_inter_dist, "\t", $avg_intra_dist, "\t", $dunn, "\n";
	}
	
	close OUT;
	#last;
	}
}

sub calc_avg_dist{
	my @gi1 = @{shift @_};
	my @gi2 = @{shift @_};
	my $dist = 0;
	my($i, $j);
	my($gi1, $gi2);
	my $n = 0;
	
	#print $#gi1+1, "\t", $#gi2+1, "\n";
	#mean-dist
	foreach $gi1 (@gi1){
		#print $gi1, "\t";
		foreach $gi2 (@gi2){
			if($gi1 != $gi2){
				#print $dist_file, "\n";
				#print $gi1, "\t", $gi2, "\t", $dists{$gi1}{$gi2}[1], "\n";
				#if($dists{$gi1}{$gi2} !~ /NA/){
					$dist += $dists{$gi1}{$gi2};
				#}
				$n++;
			}
		}
	}
	$n = 1 if $n == 0;
	
	$dist /= $n;
	
	return $dist;
}

sub calc_dunn{
	my $dist_cut = shift @_;
	my @clusts = @{shift @_};
	my $intra_dist;
	my $inter_dist;
	my $max_intra_dist = 0;
	my $min_inter_dist = 100;
	my $avg_intra_dist = 0;
	my $avg_inter_dist = 0;
	my $dunn = 0;
	my @gi1;
	my @gi2;

	
	my($i, $j);
	
	for($i = 0; $i <= $#clusts; $i++){
		$intra_dist = 0;
		@gi1 = keys %{$clusts{$dist_cut}{$clusts[$i]}};
		
		for($j = $i; $j <= $#clusts; $j++){
			$inter_dist = 0;
			@gi2 = keys %{$clusts{$dist_cut}{$clusts[$j]}};
				
				if($i == $j){
					$intra_dist = calc_avg_dist(\@gi1, \@gi2);
					$max_intra_dist = $intra_dist if $max_intra_dist < $intra_dist && $intra_dist != 0;
					$avg_intra_dist += $intra_dist/($#clusts+1);
					#print $clusts[$i], "\t", $clusts[$j], "\t", $intra_dist;
					#print "\t*", $max_intra_dist if $intra_dist == $max_intra_dist;
					#print "\n";
					#print $#gi1+1, "\t", $#gi2+1, "\n"; 
				}
				else{
					$inter_dist = calc_avg_dist(\@gi1, \@gi2);
					$min_inter_dist = $inter_dist if $min_inter_dist > $inter_dist && $inter_dist != 0;
					$avg_inter_dist += $inter_dist/(($#clusts+1)*($#clusts)/2);
					#print $clusts[$i], "\t", $clusts[$j], "\t", $inter_dist; 
					#print "\t*", $min_inter_dist, "\n" if $min_inter_dist == $inter_dist;
					#print "\n";
					#print $#gi1+1, "\t", $#gi2+1, "\n"; 
				}
		}
	}
		
	$max_intra_dist = 1 if $max_intra_dist == 0;
	$min_inter_dist = 1 if $min_inter_dist == 100;

	$dunn = $min_inter_dist/$max_intra_dist;
	
	return ($avg_inter_dist, $avg_intra_dist, $dunn);
}

sub calc_silhouette{
	my $dist_cut = shift @_;
	my @clusts = @{shift @_};
	my $num_seqs;
	my $num_clusts = $#clusts + 1;
	my($dist_intra, $dist_inter);
	my $dist;
	my $gi_ref;
	my $gi;
	
	my @gi1;
	my @gi2;
	my %sil = ();
	my $silhouette = 0;
	my $sil_width;

	my %d = (); # average intra and inter- cluster distance hash = ( cluster -> protein -> [intra_clust avg. distance, lowest average inter_clust distance, neighbouring cluster])
	my($i, $j);
	
	for($i = 0; $i <= $#clusts; $i++){
		@gi1 = keys %{$clusts{$dist_cut}{$clusts[$i]}};
		#print $dist_cut, "\t", $clusts[$i], "\t", join "\t", @gi1, "\n";
		#print $dist_cut, "\t", $clusts[$i];
		foreach $gi (@gi1){
			$d{$clusts[$i]}{$gi} = [0,0,0];
			#print "\t", $gi;
			for($j = 0; $j <= $#clusts; $j++){
				$dist = 0;
				if($j == $i){
					#calculate average intra_clust silhouette for each protein
					$d{$clusts[$i]}{$gi}[0] = calc_avg_dist([$gi], \@gi1);
				}
				else{
					#calculate average inter-clust silhouette for each protein
					@gi2 = keys %{$clusts{$dist_cut}{$clusts[$j]}};
					$dist = calc_avg_dist([$gi], \@gi2);
					
					#store the lowest average inter-cluster distance of gi (nearest-neighbour cluster)
					if($d{$clusts[$i]}{$gi}[0] == 0){
						$d{$clusts[$i]}{$gi}[1] = $dist;
						$d{$clusts[$i]}{$gi}[2] = $clusts[$j];
					}
					elsif($d{$clusts[$i]}{$gi}[0] > $dist){
						$d{$clusts[$i]}{$gi}[1] = $dist;
						$d{$clusts[$i]}{$gi}[2] = $clusts[$j];
					}
				}
			}
		}
		#print "\n";
	}

	my $clust;
	my $clust_n;
	
	#calculate the average silhouette width for the given set of clusters
	$sil_width = 0;
	
	
	foreach $clust (keys %d){
		$num_seqs = keys %{$d{$clust}};
		foreach $gi (keys %{$d{$clust}}){
			($dist_intra, $dist_inter, $clust_n) = @{$d{$clust}{$gi}};
			
			
			if($dist_intra < $dist_inter){
				$silhouette = 1 - ($dist_intra/$dist_inter);
			}
			elsif($dist_intra > $dist_inter){
				$silhouette = ($dist_inter/$dist_intra) - 1;
			}
			elsif($dist_intra == $dist_inter){
				$silhouette = 0;
			}
			$sil{$clust}{$gi} = [$silhouette, $clust_n, $dist_intra, $dist_inter];
			$sil_width += $silhouette/$num_seqs;
		}
	}
	
	#print out silhouette plots
	if(0){
		
		print "###dist_cut: $dist_cut\n";
		print "#cluster\tneighbour_cluster\tgi\tavg_intra_cluster_distance\tavg_inter_cluster_distance\tsilhouette\n";
		foreach $clust (keys %sil){
			foreach $gi (sort {$sil{$clust}{$b}[0] <=> $sil{$clust}{$a}[0]} keys %{$sil{$clust}}){
				($silhouette, $clust_n, $dist_intra, $dist_inter) = @{$sil{$clust}{$gi}};
				print join "\t", $clust, $clust_n, $gi, $dist_intra, $dist_inter, $silhouette, "\n";
			}
		}
		print "###avg_silhouette_width = ", $sil_width/$num_clusts, "\n\n";
	}

	return $sil_width/$num_clusts;
}

sub calc_simpson{
	my $dist_cut = shift @_;
	my $tot_seq = shift @_;
	my $num_clusts = shift @_;

	my $recip_simp;
	my $simp = 0;
	my @gi = ();
	my $clust_size;
	
	foreach $clust_num (keys %{$clusts{$dist_cut}}){
		@gi = keys %{$clusts{$dist_cut}{$clust_num}};
		$clust_size = $#gi + 1;
		
		$simp +=  $clust_size**2;
	}
	$simp /= $tot_seq**2;
	
	$recip_simp = $num_seqs if $num_clusts == $tot_seq;
	$recip_simp = 1/$simp if $num_clusts != $tot_seq && $simp != 0;
	
	return $recip_simp;

}

sub calc_shannon {
	my $dist_cut = shift @_;
	my $num_clusts = shift @_;
	my @clusts = @{shift @_};
	my $tax_cut = shift @_;
	
	my($intra_shannon, $inter_shannon) = (0,0);
	my $s;
	my $gi;
	my @gi1;
	my @gi2;
	my @tax;
	my $t;
	my %t;
	my $num_phyla;
	my $total_phyla;
	my $clust_pairs = $num_clusts*($num_clusts-1)/2;
	$clust_pairs = 1 if $num_clusts == 1;
	
	my($i, $j);

	
	for($i = 0; $i <= $#clusts; $i++){
		
		@gi1 = keys %{$clusts{$dist_cut}{$clusts[$i]}};
		
		for($j = $i; $j <= $#clusts; $j++){
			
			@gi2 = keys %{$clusts{$dist_cut}{$clusts[$j]}};
			%t = ();
			$total_phyla = 0;
			foreach $gi (@gi1){
				if(exists $tax_dist{$gi}){
					@tax = split /; /, $tax_dist{$gi};
					
					$t = $tax[$tax_cut] if $#tax >= $tax_cut;
					$t = "NA" if $#tax < $tax_cut;
					$t{$t}++ if exists $t{$t};
					$t{$t} = 1 if !exists $t{$t};
				}
				else{
					$t{"NA"} = 1;
				}
				$total_phyla++;
			}
			
			if($i != $j){
				foreach $gi (@gi2){
					if(exists $tax_dist{$gi}){
						@tax = split /; /, $tax_dist{$gi};
						
						$t = $tax[$tax_cut] if $#tax >= $tax_cut;
						$t = "NA" if $#tax < $tax_cut;
						$t{$t}++ if exists $t{$t};
						$t{$t} = 1 if !exists $t{$t};
					}
					else{
						$t{"NA"} = 1;
					}
					$total_phyla++;
				}	
			}
			
			foreach $t (keys %t){
				#if($t ne "NA"){
					$s = -($t{$t}/$total_phyla)*log($t{$t}/$total_phyla);
					if($i == $j){
					
						#print $clusts[$i], "\t", $t, "\t", $t{$t}, "\n";
						$intra_shannon += $s;
					}
					else{
						$inter_shannon += $s;
					}
				#}
			}
		}
	}
	
	return ($intra_shannon/$num_clusts, $inter_shannon/$clust_pairs);


}