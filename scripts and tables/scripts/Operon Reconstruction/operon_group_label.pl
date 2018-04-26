#operon_group_label.pl

use strict;
use warnings;

my $choose = shift @ARGV;
my @c = ("cellulose", "acet-cellulose", "pel", "pnag", "alginate");

while(!$choose){
	print "Please indicate which operon you would like to parse (cellulose, acet-cellulose, pel, pnag, or alginate): ";
	
	$choose = <>;
	chomp $choose;
}


my $operon_file = "./$choose\_op_final.txt";

my $outfile;

my $phylo_group_file;
my @phylo_group_files;

my @liner;
my %groups;
my %groups2;
my $seq;
my $g;
my $s;
my @s;
my $sub;
my @sub;
my %na_num;

my $cut_set;
my @cut_set = (0, 1, 2, 3);



foreach $cut_set (@cut_set){
	$outfile = "$root_dir/analysis/clust_stats/operons/$choose/$choose\_operons_phylogroups_$cut_set.txt";
	print $outfile, "\n";
	
	@phylo_group_files = <"$root_dir/analysis/clust_stats/groups/$choose/*_groups_$cut_set.txt">;
	
	%groups = ();
	%groups2 = ();
	foreach $phylo_group_file (@phylo_group_files){
		$sub = $1 if $phylo_group_file =~ /$choose\/(.*)_groups/;
		
		$na_num{$sub} = 0;
		
		open GROUP, $phylo_group_file or die;

		while(<GROUP>){
			chomp;
			@liner = split /\t/;
			$g = splice @liner, 0, 1;

			foreach $seq (@liner){
				$groups{$seq}{$sub} = $g;
				$groups2{$seq} = $g;
				
				#print $seq, "\n", $g, "\n" if $seq =~ /YP_007708745/;
			}
		}
		close GROUP;
	}	

	
	my %get_seqs;
	my $g_found;
	my $gi;
	my $op_group;
	my $group;
	my $num;

	open OUT, ">", $outfile;
	open OP, $operon_file;

	while(<OP>){
		chomp;
		@liner = split /\t/;
	
		$g_found = "";
	
		foreach $seq (@liner[4..$#liner]){
			@s = split /\|/, $seq;
			if($#s > 0){
				#print $liner[2], "\t", $seq, "\t", $#s, "\n";
				#print $s[$#s-2], "\n";
				$s = $s[1];
				$group = "";
			
				if(0){
					@sub = split /\+/, $s[0];
					$num = $#sub + 1;
					#print $s[0], "\n" if $s[0] =~ /\+/;
				
					foreach $sub (@sub){
				
						$num--;
						if(exists $groups{$s}{$sub}){
					
							$g = $groups{$s}{$sub};
							$get_seqs{$s} = [$g, @liner[2..3]];
					
							$group .= $g . "($s)";
							$group .= "+" if $num;

						}
						else{
							$na_num{$sub}++;
						
							$g = $sub . "_NA" . $na_num{$sub};
							$get_seqs{$s} = [$g, @liner[2..3]];
					
							$group .= $g . "($s)";
							$group .= "+" if $num;
						}
					}
				}
				else{
					if(exists $groups2{$s}){
					
						$g = $groups2{$s};
						$get_seqs{$s} = [$g, @liner[2..3]];
					
						$group .= $g . "($s)";

					}
					else{
						$na_num{$s[0]}++;
						$g = $s[0] . "_NA" . $na_num{$s[0]};
						$get_seqs{$s} = [$g, @liner[2..3]];
					
						$group .= $g . "($s)";
					}
				}
					
				$g_found .= "|" . $group if $g_found ne ""; 
				$g_found = $group if $g_found eq "";
				#print $g_found, "\n" if $s[0] =~ /\+/;
			}
		}
		if($g_found){# !~ /NA/){
			print OUT $g_found, "\t", join "\t", @liner;
			print OUT "\n";
		}
		elsif(0){
			print $g_found, "\t", join "\t", @liner, "\n";
		}
	}

	close OUT;
	close OP;

}