use ExtUtils::MakeMaker;

WriteMakefile(
	NAME          => 'Microarray::GEO::SOFT',
	VERSION_FROM  => 'lib/Microarray/GEO/SOFT.pm',
	PREREQ_PM     => {
		'List::Vectorize' => "",
		'Microarray::ExprSet' => "",
		'File::Basename' => "",
		'LWP::UserAgent' => "",
		'Thread' => "",
		'threads::shared' => "",
		'Time::HiRes' => "",
		'Test::More' => "",
	},
	AUTHOR        => 'Zuguang Gu <jokergoo@gmail.com>',
	ABSTRACT      => 'Reading microarray data in SOFT format from GEO database',
	LICENSE       => 'perl',
);
