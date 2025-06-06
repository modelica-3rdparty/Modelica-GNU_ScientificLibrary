within GNU_ScientificLibrary.Functions.integration;

impure function integrand_setup "Integrand Setup for GSL integration routines"
  // contact: tommy.burch@physik.uni-r.de
  extends Modelica.Icons.Function;
  input String C_integrand = "exp(-x*x)";
  output Integer rcode[2];
  protected String it_file,mi_file,mi_hfile,ar_file,incl_path;
  protected String C_line = "f=" + C_integrand + ";return(f);}";

algorithm
// need system commands to construct the integrand C function, compile the object file,
//   and build the archive file
//  cwd := Modelica.Utilities.System.getWorkDirectory();
  it_file := Modelica.Utilities.Files.loadResource("modelica://GNU_ScientificLibrary/Resources/Source/itop.c");
  mi_file := Modelica.Utilities.Files.loadResource("modelica://GNU_ScientificLibrary/Resources/Source/integration_mi.c");
  mi_hfile := Modelica.Utilities.Files.loadResource("modelica://GNU_ScientificLibrary/Resources/Include/integration_mi.h");
  incl_path := Modelica.Utilities.Files.fullPathName("modelica://GNU_ScientificLibrary/Resources/Include");
  Modelica.Utilities.Files.copy(it_file,"integrand.c",replace=true);
  Modelica.Utilities.Files.copy(mi_file,"integration_mi.c",replace=true);
  Modelica.Utilities.Files.copy(mi_hfile,"integration_mi.h",replace=true);
  Modelica.Utilities.Streams.print(C_line,"integrand.c");
  rcode[1] := Modelica.Utilities.System.command("gcc -fPIC -c integration_mi.c integrand.c -I" + incl_path);
  rcode[2] := Modelica.Utilities.System.command("gcc -shared integration_mi.o integrand.o -lm -llibgsl -o libgsl_integration_MI.so");
// pointing out other possible LibraryDirectory (e.g., in qag) seems not to work:
//   copying archive file to default directory...
  ar_file := Modelica.Utilities.Files.loadResource("modelica://GNU_ScientificLibrary/Resources/Library/linux64/liblibgsl_integration_MI.so");
  Modelica.Utilities.Files.copy("libgsl_integration_MI.so",ar_file,replace=true);
  ar_file := Modelica.Utilities.Files.loadResource("modelica://GNU_ScientificLibrary/Resources/Library/win32/libgsl_integration_MI.dll");
  Modelica.Utilities.Files.copy("libgsl_integration_MI.so",ar_file,replace=true);
  ar_file := Modelica.Utilities.Files.loadResource("modelica://GNU_ScientificLibrary/Resources/Library/win64/libgsl_integration_MI.dll");
  Modelica.Utilities.Files.copy("libgsl_integration_MI.so",ar_file,replace=true);
  annotation(
    Documentation(info = "<html><head></head><body>Integrand-setup \"function\" for using GSL integration routines within a larger Modelica model.<div><br></div><div>Wherever possible, Modelica.Utilities routines are used to handle file / system commands. The current exception is the use of 'gcc' to compile the libgsl_integration_MI.so / libgsl_integration_MI.dll library. This is standard on any Linux system and is available for any Windows system via 'mingw'.</div><div><br></div><div><div>Alternatively, one could edit integrand.c directly and then compile \"by hand\" in a command window:&nbsp;</div><div><span class=\"Apple-tab-span\" style=\"white-space:pre\">	</span>gcc -fPIC -c integrand.c integration_mi.c</div><div><span class=\"Apple-tab-span\" style=\"white-space:pre\">	</span>gcc -shared integration.o integration_mi.o -o libgsl_integration_MI.so &nbsp;[.dll for windows]</div></div></body></html>"));
end integrand_setup;