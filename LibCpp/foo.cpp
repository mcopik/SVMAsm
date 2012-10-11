char string[] = "Ala ma kota";

extern "C" int foo(char * x) {
	*x = string[0];
	return 5;
}
