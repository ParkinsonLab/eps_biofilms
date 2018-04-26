#groups_relabel.pl

use strict;
use warnings;

my $sys_chosen = shift @ARGV;
my @op_files = <"./operons/$sys_chosen\_*">; #Change to the original groups file??
my @dist_files = <"./dists/$sys_chosen/*_dists.txt">;

my $outfile = "./group_relabel/$sys_chosen\_groups_relabel.txt";

my($op_file, $dist_file);

my $sub;
my @liner;
my @g;
my($i, $j);
my($grp1, $grp2, $gi1, $gi2);
my %grp_comb;
my %grp_count;
my %gi_grp;

my %dist;
my $dist;
my $n;
my %grp_order;
my @g_order;

my @g2;
my $match;
my %grp_match;
my %match_max;
my %grp_tmp;
my %merge_grp;

my $grp_new;
my $grp_num;
my $grp_suffix;
my $na_nosuff;
my $na_num;
my $g_num;

my %grp_assign;


print $outfile, "\n";

open OUT, ">", $outfile;

print OUT "#old_group\tnew_group\n";


foreach $op_file (@op_files){
	
	open OP, $op_file;
	
	while(<OP>){
		chomp;
		if($. != 1){
			@liner = split /\t/;
			
			@g = split /\|/, $liner[0];
			
			for($i = 0; $i <= $#g; $i++){
				($grp1, $gi1) = ($1, $2) if $g[$i] =~ /^(.*)\((\d*)\)/;
				$sub = $1 if $grp1 =~ /(.*)_/;

				for($j = 0; $j <= $#g; $j++){
					($grp2, $gi2) = ($1, $2) if $g[$j] =~ /^(.*)\((\d*)\)/;
					$grp_comb{$sub}{$grp1}{$grp2} = 1 if $grp1 ne $grp2;
				}
				
				$grp_count{$sub}{$grp1}{$gi1} = 1;
				$gi_grp{$sub}{$gi1} = $grp1;
			}
						
		}
	
	}
	close OP;
	
	#order groups by phylogenetic distance to the largest sequence cluster
	foreach $dist_file (@dist_files){
		$sub = $1 if $dist_file =~ /$sys_chosen\/(.*)_dists/;
	
		open DIST, $dist_file;
		
		while(<DIST>){
			chomp;
			@liner = split /\t/;
			($gi1, $gi2, $dist) = @liner[0..2];
			
			$dist{$sub}{$gi1}{$gi2} = $dist;
			$dist{$sub}{$gi2}{$gi1} = $dist;

		}
		
		close DIST;
		
		@g = sort {keys %{$grp_count{$sub}{$b}} <=>	keys %{$grp_count{$sub}{$a}}} keys %{$grp_count{$sub}};
		$grp1 = shift @g;
		$grp_order{$sub}{$grp1} = 0;
		
		for($i = 0; $i <= $#g; $i++){
			$grp2 = $g[$i];
			
			$dist = 0;
			$n = 0;
			
			foreach $gi1 (keys %{$grp_count{$sub}{$grp1}}){
				foreach $gi2 (keys %{$grp_count{$sub}{$grp2}}){
					if(exists $dist{$sub}{$gi1}{$gi2}){
						$dist += $dist{$sub}{$gi1}{$gi2};
						$n++;
					}
				}
			}
			
			if($n > 0){
				$grp_order{$sub}{$grp2} = $dist/$n;
			}
			else{
				#arbitrary token indicates if an operon-associated locus is an un-clustered hit found by iterative hmm searches
				$grp_order{$sub}{$grp2} = 1000;
			}
		}
	
	}
	
	#exit;
	
	%grp_match = ();
	
	foreach $sub (sort {$a cmp $b} keys %grp_comb){
		@g = sort {$a cmp $b} keys %{$grp_comb{$sub}};

		for($i = 0; $i <= $#g; $i++){
			$grp1 = $g[$i];
			$match_max{$sub}{$grp1} = 0;
			@g2 = sort {$a cmp $b} keys %{$grp_comb{$sub}{$grp1}};
			
			for($j = 0; $j <= $#g; $j++){
				$grp2 = $g[$j];
				
				if($grp1 ne $grp2){
				
					$match = grp_match(\@g2, $grp2, $sub);
					
					#print join "\t", $grp1, @g2, "\n" if $grp1 =~ /Alg44_G3|(Alg44_G3|Alg44_G9|Alg44_G13)/ && $grp2 =~ /Alg44_G3|(Alg44_G3|Alg44_G9|Alg44_G13)/ && $sub eq "Alg44";
				
					if($match > 1){
						$grp_match{$sub}{$grp1}{$grp2} = $match;
					
						if($match > $match_max{$sub}{$grp1}){
							$match_max{$sub}{$grp1} = $match;
						}
					}
				}
			}
		}
		
		#print $match_max{"Alg44_G3"}, "\n";
		
		#create higher-order sequence cluster groupings for a given EPS locus group based on the 
		#maximum number of unique co-associated EPS sequence clusters shared at other operon loci.
		%grp_tmp = ();
		$n = 1;
		foreach $grp1 (sort {$a cmp $b} keys %{$grp_match{$sub}}){

			foreach $grp2 (keys %{$grp_match{$sub}{$grp1}}){
				$match = $grp_match{$sub}{$grp1}{$grp2};
				
				if($match == $match_max{$sub}{$grp1} && $match >= $match_max{$sub}{$grp2}){
					
					if(!exists $grp_tmp{$grp1} && !exists $grp_tmp{$grp2}){
						$grp_tmp{$grp1} = $n;
						$grp_tmp{$grp2} = $n;
						$n++;	
					}
					elsif(exists $grp_tmp{$grp1} && !exists $grp_tmp{$grp2}){
						$grp_tmp{$grp2} = $grp_tmp{$grp1};
					
					}
					elsif(!exists $grp_tmp{$grp1} && exists $grp_tmp{$grp2}){
						$grp_tmp{$grp1} = $grp_tmp{$grp2};
					}
						
					#print $grp1, "\t", $grp2, "\t", $match, "\n";
				}
			}	
		}
		
		
		
		
		foreach $grp1 (sort {$grp_tmp{$a} <=> $grp_tmp{$b}} keys %grp_tmp){
			$n = $grp_tmp{$grp1};
			$merge_grp{$sub}{$n}{$grp1} = 1;
			#print $grp1, "\t", $grp_tmp{$grp1}, "\t", $grp_order{$sub}{$grp1}, "\n";
		}
		#print "####\n";
		
		
		
		#For each EPS subunit/locus:
		#Relabel sequence groups according to mean phylogenetic distance to the sequence group with the greatest number of sequences; 
		#Also merge groups into a higher-order grouping label if they are also share a maximum number of unique sequence clusters at other loci.
		
		
		@g_order = sort {$grp_order{$sub}{$a} <=> $grp_order{$sub}{$b}} keys %{$grp_order{$sub}};
		
		$grp_num = 1;
		$na_num = 1;
		$g_num = 1;
		$na_nosuff = 1;
				
		for($i = 0; $i <= $#g_order; $i++){
			$grp1 = $g_order[$i];
			
			#print $grp1, "\t", $grp_order{$sub}{$grp1}, "\n";
			
			$n = keys %{$merge_grp{$sub}};
			
			if(exists $grp_tmp{$grp1} && !exists $grp_assign{$sub}{$grp1} && $n > 1){
				$n = $grp_tmp{$grp1};
				
				@g = keys %{$merge_grp{$sub}{$n}};
				$na_num = 1;
				$g_num = 1;
					
				foreach $grp2 (sort {$grp_order{$sub}{$a} <=> $grp_order{$sub}{$b}} @g){
					if($grp_order{$sub}{$grp2} == 1000 || $grp2 =~ /NA/){
						$grp_suffix = "NA" . $na_num;
						$na_num++;
					}
					else{
						$grp_suffix = $g_num;
						$g_num++;
					}
					$grp_new = $sub . "_G" . $grp_num . "-" . $grp_suffix;
						
					$grp_assign{$sub}{$grp2} = $grp_new;
						
				}
				$grp_num++;
				
			}
			elsif($n == 1 || !exists $grp_assign{$sub}{$grp1}){
				if($grp_order{$sub}{$grp1} == 1000 || $grp1 =~ /NA/){
						$grp_suffix = "_NA" . $na_nosuff;
						$na_nosuff++;
				}else{
					$grp_suffix = "_G" . $grp_num;
					$grp_num++;
				}
				
				$grp_new = $sub . $grp_suffix;
				$grp_assign{$sub}{$grp1} = $grp_new;
			}
			
			#print $grp1, "\t", $grp_order{$sub}{$grp1}, "\n";
		}
		
		
		###Finally, output the relabelled groups
		
		foreach $grp1 (sort {$grp_order{$sub}{$a} <=> $grp_order{$sub}{$b}} keys %{$grp_assign{$sub}}){
			print OUT $grp1, "\t", $grp_assign{$sub}{$grp1}, "\n";
			#print $grp1, "\t", $grp_order{$sub}{$grp1}, "\t", $grp_assign{$sub}{$grp1}, "\n";
		}

		#last;
	}
	
}

close OUT;

sub grp_match{
	my @g2 = @{$_[0]};
	my($grp2, $sub) = @_[1..2];
	my $g;
	my $grp_match = 0;

	
	foreach $g (@g2){
		$grp_match++ if exists $grp_comb{$sub}{$grp2}{$g};
	}
	
	return $grp_match;
	
}