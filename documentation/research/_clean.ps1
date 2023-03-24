
$types = ".aux", ".bbl", ".log", ".toc", ".blg"

foreach ($type in $types)
{
	Remove-Item *$type
}
