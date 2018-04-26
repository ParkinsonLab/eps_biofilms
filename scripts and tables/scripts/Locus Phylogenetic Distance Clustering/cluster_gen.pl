#cluster_gen.pl

use strict;
use warnings;

my @infiles = <"./phylogenetic_distances/*/*_dists.txt">;
my $infile;
my $outfile;
my $outdir = "./tmp/";

my $dist_cut = 1;
my $max_dist;
my $i;

my($gi, $gi1, $gi2, $dist, $dist1, $dist2);

my($clust_num, $clust_num1, $clust_num2);
my $clust_last;
my %clust_num;
my %clusts;
my %clusts_tmp;
my %dists;

my @gi;

my $line;
my @in;

foreach $infile (@infiles){
	
	#if($infile =~ /Wss/){
	$outfile = $infile;
	$outfile =~ s/^.*\///;
	$outfile = $outdir . $outfile;
	$outfile =~ s/_dists/_clusters/;

	
	$max_dist = 0;
	@in = ();
	
	open IN, $infile;
	
	while(<IN>){
		chomp;
		($gi1, $gi2, $dist) = split /\t/;
		$dists{$gi1}{$gi2} = $dist;
		$dists{$gi2}{$gi1} = $dist;
		$max_dist = $dist if $dist > $max_dist;
		push @in, $_;
	}
	
	close IN;
	
	#print join "\n", @in, "\n";
	%clusts = ();
	%clusts_tmp = ();
	%clust_num = ();
	for($dist_cut = 0; $dist_cut <= $max_dist+0.5; $dist_cut += 0.01){
	#$dist_cut = 2.17;
	#foreach $dist_cut (2.17){

		($clust_num1, $clust_num2) = (0, 0);
		$clust_last = 0;
		
		foreach $line (@in){
			#print $line, "\n";
			($gi1, $gi2, $dist) = split /\t/, $line;

		
			if($dist <= $dist_cut){
				($clust_num1, $clust_num2) = (0, 0);
				($dist1, $dist2) = (1000, 1000);
			
				($clust_num1, $dist1) = @{$clust_num{$dist_cut}{$gi1}} if exists $clust_num{$dist_cut}{$gi1};
				($clust_num2, $dist2) = @{$clust_num{$dist_cut}{$gi2}} if exists $clust_num{$dist_cut}{$gi2};
			
				if($clust_num1 != $clust_num2){
					if($clust_num1 == 0){
						$clust_num{$dist_cut}{$gi1} = [$clust_num2, $dist];
						$clusts_tmp{$dist_cut}{$clust_num2}{$gi1} = 1;
				
					}
					elsif($clust_num2 == 0){
						$clust_num{$dist_cut}{$gi2} = [$clust_num1, $dist];
						$clusts_tmp{$dist_cut}{$clust_num1}{$gi2} = 1;
					}
					else{
						if($dist < $dist2){
							$clust_num{$dist_cut}{$gi2} = [$clust_num1, $dist];
							$clusts_tmp{$dist_cut}{$clust_num1}{$gi2} = 1;
							delete $clusts_tmp{$dist_cut}{$clust_num2}{$gi2};
					
						}
						elsif($dist < $dist1){
							$clust_num{$dist_cut}{$gi1} = [$clust_num2, $dist];
							$clusts_tmp{$dist_cut}{$clust_num2}{$gi1} = 1;
							delete $clusts_tmp{$dist_cut}{$clust_num1}{$gi1};	
						}
					}
				}
				elsif($clust_num1 == 0 && $clust_num2 == 0){
					$clust_last++;
					$clust_num = $clust_last;
					$clust_num{$dist_cut}{$gi1} = [$clust_num, $dist];
					$clust_num{$dist_cut}{$gi2} = [$clust_num, $dist];
				
					$clusts_tmp{$dist_cut}{$clust_num}{$gi1} = 1;
					$clusts_tmp{$dist_cut}{$clust_num}{$gi2} = 1;
				}
			}
			else{
				foreach $gi ($gi1, $gi2){
					if(!exists $clust_num{$dist_cut}{$gi}){
						$clust_last++;
						$clust_num{$dist_cut}{$gi} = [$clust_last, 1000];
						$clusts_tmp{$dist_cut}{$clust_last}{$gi} = 1;
					}
				}
			}
		}
		#print $dist_cut, "\n";
		
		$clust_last = 0;
		
		foreach $clust_num (sort {keys %{$clusts_tmp{$dist_cut}{$b}} <=> keys %{$clusts_tmp{$dist_cut}{$a}}} keys %{$clusts_tmp{$dist_cut}}){
			@gi = keys %{$clusts_tmp{$dist_cut}{$clust_num}};
			
			if($#gi != -1){
				$clust_last++;
				foreach $gi (keys %{$clusts_tmp{$dist_cut}{$clust_num}}){
					$clusts{$dist_cut}{$clust_last}{$gi} = 1;
					#print $gi, "\t", $clust_last, "\t", $clust_num, "\t", $clust_num{$dist_cut}{$gi}[0], "\n";
					$clust_num{$dist_cut}{$gi}[0] = $clust_last;
				}
			}
		}
		
		last if $clust_last == 1;
	}
	
	
	@gi = keys %{$clust_num{0}};
	my @dists = sort {$a <=> $b} keys %clust_num;
	#print $max_dist;
	
	print $outfile, "\n";
	
	if(1){
		open OUT, ">", $outfile;
		print OUT "\t", join "\t", @dists, "\n";
	
		foreach $gi (@gi){
			print OUT $gi;
			foreach $dist_cut (@dists){
				print OUT "\t", $clust_num{$dist_cut}{$gi}[0];
			}
			print OUT"\n";
		}
		
		close OUT;
	}
	elsif(0){
		my $j;
		print $dist_cut, "\n";
		foreach $clust_num (sort {keys %{$clusts{$dist_cut}{$b}} <=> keys %{$clusts{$dist_cut}{$a}}} keys %{$clusts{$dist_cut}}){
			
			
			@gi = keys %{$clusts{$dist_cut}{$clust_num}};
			for($i = 0; $i < $#gi; $i++){
				for($j = $i+1; $j <= $#gi; $j++){
					#print $clust_num, "\t", $gi[$i], "\t", $gi[$j], "\t", $dists{$gi[$i]}{$gi[$j]}, "\n";
				}
			}
			#print "clust_", $clust_num, "\t", $#gi+1, "\t";
			#print join "\t", @gi, "\n";
		}
	}
	
	#print $clust_last, "\n";
	#}
	last;
}
