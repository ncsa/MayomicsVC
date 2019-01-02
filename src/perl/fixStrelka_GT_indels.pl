my @nfields;
my @samples;
while(<>)       {
        if (/^\#\#/) {
                print $_;
                next;
        }
        if (/^\#/) {
                print "##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">\n";
                chomp;
                my $line = $_;
                $line =~ s/\#//;
                @nfields = split (/\t/, $line) ;
                @samples = @nfields[9..$#nfields];
                print "#$line\n";
                next;
        }
        chomp;
        my %fields;
        my @values = split/\t/;
        
        foreach my $i (@nfields) {
                        $fields{$i} = shift @values;
        }
        my @names_format = split (/:/, $fields{'FORMAT'});
        push(@names_format, 'GT');
	my @values_info = split (/;/, $fields{'INFO'});
	### for indels strelka a different annotation in SGT field
	#ref,het,hom,conflict
	my ($NormalAllele,$TumorAllele);		
	foreach my $i (@values_info)	{
		if ($i =~ /^SGT=/)	{
			my ($sgt,$gtvalue)=split('=',$i);
			($NormalAllele,$TumorAllele)=split('->',$gtvalue);
			if($NormalAllele eq 'ref')	{
				$NormalAllele="0/0";
			} elsif($NormalAllele eq 'het')	{
				$NormalAllele="0/1";
			} elsif($NormalAllele eq 'hom')	{
				$NormalAllele="1/1";
			}else	{$NormalAllele="./.";}	

			if($TumorAllele eq 'ref')      {
                                $TumorAllele="0/0";
                        } elsif($TumorAllele eq 'het') {
                                $TumorAllele="0/1";
                        } elsif($TumorAllele eq 'hom') {
                                $TumorAllele="1/1";
                        }else   {$TumorAllele="./.";}
		}
	}
	my (%sample_values);
	foreach my $sample (@samples) {
		my @values;
		my @s_values = split (/:/, $fields{$sample});
		foreach my $i (@names_format) {
			my $cvalue = shift (@s_values);
			$sample_values{$sample}->{$i}=$cvalue;
			if ($i eq 'GT')	{
				if ($sample eq 'NORMAL')	{
					push(@values, $NormalAllele);
				}else { 
					push(@values,$TumorAllele);
				}
			}else	{
				push(@values, $sample_values{$sample}->{$i});
			}
		}
	@values=reverse @values;
	$fields{$sample} = join (':', @values);

	}
@names_format= reverse @names_format;
	
my @line_fields;
        $fields{'FORMAT'} = join (':', @names_format);
        foreach my $field (@nfields) {
        push (@line_fields, $fields{$field});
    }
    print join ("\t",@line_fields) . "\n";
}





