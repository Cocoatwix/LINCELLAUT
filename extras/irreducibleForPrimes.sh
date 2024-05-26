primes=(2 3 5 7 11 13 17 19 23 29 31 37 41 43 47)
newPrimes=(2)
sizeOfPrimes=15
isPrime=0
upperBound=70
echo "Getting prime numbers up to $upperBound..."
for ((i = 50; i < $upperBound; i += 1)); do
	isPrime=1
	for p in ${primes[@]}; do
		if [ $(($i % $p)) -eq 0 ]; then
			isPrime=0
			break
		fi
	done
	if [ $isPrime -eq 1 ]; then
		primes+=($i)
		newPrimes+=($i)
		sizeOfPrimes=$(($sizeOfPrimes+1))
	fi
done
echo "Prime numbers found: ${newPrimes[*]}. Running irreducible with each prime number as the modulus..."
for n in ${newPrimes[@]}; do
    ./irreducible $n
done
echo "Done."