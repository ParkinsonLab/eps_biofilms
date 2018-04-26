#neighbour_network.pl - generate a genome proximity network (multi-partite) of putative EPS subunit hits
#Question: which hits/sequence families are frequently found in proximity to each other in bacterial genomes?

#Calculate the average genome proximity (in KBs) between genes encoding putative EPS subunit hits 
#that belong to distinct sequence groups (95% needleman-wunsch pairwise alignment sequence identity).
#Do sequence groups for distinct EPS subunit contain sequences from the same set of genomes?


use strict;
use warnings;

my $root_dir = "E:/eps_latest";

#my @operon_files = <"$root_dir/reconstruct_new/operons_final/*_operons_phylogroups.txt">;
my @operon_files = <"$root_dir/analysis/op_fill/operons_filled_new/*_operon_phylogroups_withfill.txt">;

my $op_ref_file = "$root_dir/data/ref_op_orders.txt";
my $infile;
my $outfile;

my @c = ("cellulose", "acet-cellulose", "alginate", "pel", "pnag");
#my $chosen = $c[2];

my $sys;
my @liner;
my %op_order;

my $chosen = shift @ARGV;

while(!$chosen){
	print "Please indicate which operon you would like to parse (cellulose, acet-cellulose, pel, pnag, or alginate): ";
	
	$chosen = <>;
	chomp $chosen;
}


my $target;
my $subunit;

my $group_num;
my $gi;
my $size;
my %seq_group;
my %avg_size;

my($sub1, $sub2);
my($s1, $e1);
my($s2, $e2);
my $dist;
my($gi1, $gi2);
my($g1, $g2);

my($sp, $tax, $sp_id);
my($loc1, $loc2);
my @t;
my $t;

my %op;
my %sp;

my %group_tally;

my %total;

my %not;

my %phyla;

my($i, $j);
my @sub;
my @gi;

my %test;

my($phylo, $phylo1, $phylo2);
my @phylo;
my %phylo_group;

my $op_file =  "$root_dir/reconstruct_new/operons_final/$chosen\_operons_phylogroups.txt";
$op_file = "$root_dir/analysis/op_fill/operons_filled_new/$chosen\_operon_phylogroups_withfill.txt";

open OPERONS, $op_file;

while(<OPERONS>){
	chomp;
	if(!/^#/){
	
		@liner = split /\t/;
	
		($phylo, $sp, $tax, $sp_id) = splice @liner, 0, 4;
		
		@t = split /; /, $tax;
		$sp{$sp_id} = [$sp, $tax];
		
		$t = $t[1] if $#t <= 4;
		$t = $t[2] if $#t > 4; #2 - phyla; 5 - genus
		
		if($#t >= 2){
			$t = $t[2] if $t[1] =~ /Proteobacteria|Firmicutes/;
			$t = "Other" if $t[1] !~ /Proteobacteria|Firmicutes/;
		}
		
		$t = "NA" if $#t < 0 || !$t;
		$t =~ s/\.//;
		
		$phyla{$t}++ if exists $phyla{$t};
		$phyla{$t} = 1 if !exists $phyla{$t};
		
	
		@phylo = split /\|/, $phylo;
		@sub = @liner[1..$#liner];
		
		
		for($i = 0; $i <= $#sub; $i++){
			($sub1, $gi1, $loc1) = split /\|/, $sub[$i];
			$phylo1 = $1 if $phylo[$i] =~ /(.*)\(/;
			
			$op{$sp_id}{$gi1} = [$phylo1, $loc1];
		
			$group_tally{$phylo1}{$t}{$sp_id} = 1;
			$total{$phylo1}{$sp_id} = 1;
		}
	}
}

close OPERONS;


foreach $sp_id (keys %op){

	@gi = sort {$a cmp $b} keys %{$op{$sp_id}};
	
	for($i = 0; $i <= $#gi; $i++){
		
		($phylo1, $loc1) = @{$op{$sp_id}{$gi[$i]}};
		($s1, $e1) = parse_loc($loc1);

		for($j = $i+1; $j <= $#gi; $j++){
			($phylo2, $loc2) = @{$op{$sp_id}{$gi[$j]}};
			($s2, $e2) = parse_loc($loc2);
				
			$dist = 0;
			if($s2 >= $e1){
				$dist = $s2 - $e1;
			}
			elsif($s1 >= $e2){
				$dist = $s1 - $e2;
			}
			else{
				$dist = $s2 - $e1;
			}
			
			if($phylo1 =~ /BcsA_G12|BcsB_G6|BcsZ_G9|BcsC_(22|21|12)/ && $phylo2 =~ /BcsA_G12|BcsB_G6|BcsZ_G9|BcsC_(22|21|12)/){
				#print $phylo1, "\t", $phylo2, "\t", $gi[$i], "\t", $gi[$j], "\t", $loc1, "\t", $loc2, "\t", $dist/1000, "\n";
			}
				
			if(!exists $phylo_group{$phylo1}{$phylo2} && !exists $phylo_group{$phylo2}{$phylo1}){
				$phylo_group{$phylo1}{$phylo2} = [1, $dist];
				#$phylo_group{$phylo2}{$phylo1} = [1, $dist];
			}
			elsif(exists $phylo_group{$phylo1}{$phylo2}){
				$phylo_group{$phylo1}{$phylo2}[0]++;
				$phylo_group{$phylo1}{$phylo2}[1] += $dist;

			}
			elsif(exists $phylo_group{$phylo2}{$phylo1}){
				$phylo_group{$phylo2}{$phylo1}[0]++;
				$phylo_group{$phylo2}{$phylo1}[1] += $dist;
			}
		}
	}
}

#exit;
my $num;
#order according to the reference operon?
#@sub = @{$op_order{$chosen}};

#Print out Genome Proximity Network?
if(1){


$outfile = "$root_dir/neighbour_network/latest/$chosen\_neighbour_network_phylo.txt";
$outfile = "$root_dir/neighbour_network/latest/withfill/$chosen\_neighbour_network_withfill_phylo.txt";

print $outfile, "\n";

open OUT, ">", $outfile;

print OUT "#Sequence_group1\tSequence_group2\tnum_occurences\tavg_genomic_dist(kb)\n";

#print out operon group associations, based around the presence of a Polysaccharide Synthase Subunit

foreach $phylo1 (keys %phylo_group){

	foreach $phylo2 (keys %{$phylo_group{$phylo1}}){
		#print $sub1, "\t", $sub2, "\n";

				
		($num, $dist) = @{$phylo_group{$phylo1}{$phylo2}};
				
		print OUT $phylo1, "\t", $phylo2, "\t", $num, "\t", ($dist/$num)/1000, "\n";
	}
}

close OUT;
}

#print out the node mapping phylogenetic distributions of Sequence Homology Groups?

if(1){

my $phy;
my @phy = sort {$phyla{$b} <=> $phyla{$a}} keys %phyla;

my $group_size;

			
$outfile = "$root_dir/neighbour_network/latest/$chosen\_network_map_phylo.txt";
$outfile = "$root_dir/neighbour_network/latest/withfill/$chosen\_network_withfill_map_phylo.txt";

print $outfile, "\n";

open OUT, ">", $outfile;
									
print OUT "#sequence_group\tsubunit_hit\ttotal_hits\t";

print OUT join "\t", @phy; 
print OUT "\n";

my $group;
my $num_seq;
my $tot;

foreach $group (sort {$a cmp $b} keys %group_tally){
	($subunit, $g1) = split /\_/, $group;
	$tot = keys %{$total{$group}};
	
	print OUT $group, "\t", $subunit, "\t", $tot;

	foreach $phy (@phy){
			
			
		if(exists $group_tally{$group}{$phy}){
			$num_seq = keys %{$group_tally{$group}{$phy}};
		}
		else{
			$num_seq = 0;
		}

		print OUT "\t", $num_seq;
	}
	print OUT "\n";
	
}

close OUT;
}

sub parse_loc{
	my $loc = shift @_;
	my $l = $loc;
	my $s;
	my $e;
	
	$l =~ s/complement//;
	$l =~ s/\(|\)//g;
	
	$l =~ s/join//;
	
	if($loc =~ /complement/){
		($s, $e) = split /\.+/, $l if $loc !~ /join/;
		
		if($loc =~ /join/){
			($s, $e) = split /,/, $l;
			$s =~ s/.*\.+//;
			$e =~ s/\.+.*//;
		}
	}
	else{
		($s, $e) = split /\.+/, $l if $loc !~ /join/;
		
		if($loc =~ /join/){
			($s, $e) = split /,/, $l;
			$s =~ s/\.+.*//;
			$e =~ s/.*\.+//;
		}
		
	}
	
	return ($s, $e);
	
}
