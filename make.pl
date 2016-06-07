use strict;

my %TOTAL_ENERGY=(
	"INTERNAL" => 0,
	#"TOTAL" => 1,
);

my %DERIVATIVE=(
	"SPH" => 0,
	"MBG" => 1,
	#"LE" => 2,
	#"MPS" => 3,
);

my %AV_LIMITTER=(
	"NO"   => 0,
	"B95"  => 1,
	"CD10" => 2,
);

my %AV_TYPE=(
	#"MG83" => 0,
	"M97" => 1,
	"vNR" => 2,
	#"LM85" => 3,
	#"NO" => 4,
);

my %TIME_VAR_AV=(
	"NO"  => 0,
	"MM97" => 1,
	"CD10" => 2,
	#"RH12" => 3,
	#"R15" => 4,
);

my $test;

foreach my $key_eng(keys %TOTAL_ENERGY){
	foreach my $key_drv(keys %DERIVATIVE){
		foreach my $key_av(keys %AV_TYPE){
			foreach my $key_shr(keys %AV_LIMITTER){
				foreach my $key_blk(keys %TIME_VAR_AV){
					my $name = "${key_eng}_${key_drv}_${key_av}_${key_shr}_${key_blk}";
					$test = $test . $name . " ";
					my $out;
					$out .= qq|$name: | . q|$(CPPSRCS) $(CPPHDRS)| . qq|\n|;
					$out .= "\t" . '@$(CPPC) $(FLAGS) $(WARNINGS) $(CPPSRCS) -o ${@}.out $(LIBS)';
					$out .= qq|-D TOTAL_ENERGY=$TOTAL_ENERGY{"$key_eng"} -D DERIVATIVE=$DERIVATIVE{"$key_drv"} -D AV_LIMITTER=$AV_LIMITTER{"$key_shr"} -D AV_TYPE=$AV_TYPE{"$key_av"} -D TIME_VAR_AV=$TIME_VAR_AV{"$key_blk"}\n|;
					$out .= "\t-mkdir $name\n";
					$out .= "\t-mkdir $name/result\n";
					$out .= "\tmv $name.out ./$name\n";
					print $out . "\n";
				}
			}
		}
	}
}

$test = "test: " . $test . "\n";

print $test;
