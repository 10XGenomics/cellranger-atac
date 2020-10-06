#ifndef _cmd_args_h
#define _cmd_args_h

template <class T>
void ExitWithMsg(std::string s, T var){
    std::cout << s << " " << var << std::endl;
    exit(EXIT_FAILURE);
}
template void ExitWithMsg<int>(std::string s, int var);
template void ExitWithMsg<double>(std::string s, double var);
template void ExitWithMsg<std::string>(std::string s, std::string var);

typedef long int int64_t;

bool directory_exists(const std::string& path)
{
    struct stat info;

    if(stat( path.c_str(), &info ) != 0)
        return false;
    else if(info.st_mode & S_IFDIR)
        return true;
    else
        return false;
}

bool file_exists(const std::string& name) {
  struct stat buffer;
  return (stat (name.c_str(), &buffer) == 0);
}

struct command_args{
    std::string MATRIX_FILE;
    std::string OUT;
    int TOPICS;
    int ITER;
    double TOL;
    int THREADS;

    void print() const{
        std::cout << "MATRIX_FILE: " << MATRIX_FILE<< std::endl;
        std::cout << "OUT: " << OUT << std::endl;
        std::cout << "TOPICS: " << TOPICS << std::endl;
        std::cout << "ITER: " << ITER << std::endl;
        std::cout << "TOL: " << TOL << std::endl;
        std::cout << "THREADS: " << THREADS << std::endl;
    }
};

typedef std::map<std::string, docopt::value> argmap;
struct itag_{};
struct dtag_{};
struct stag_{};
struct btag_{};
template <typename T> struct dc_tag;
template < > struct dc_tag<int> {typedef itag_ type;};
template < > struct dc_tag<double> {typedef dtag_ type;};
template < > struct dc_tag<std::string> {typedef stag_ type;};
template < > struct dc_tag<bool> {typedef btag_ type;};

void process_arg(int& var, argmap& args, std::string key, itag_){
     docopt::value& V = args.at(key);
     if (!V) ExitWithMsg("missing arg values for", key);
     try{ var = V.asLong(); } catch (...) { ExitWithMsg("Expected integer for", key);}
}
void process_arg(double& var, argmap& args, std::string key, dtag_){
     docopt::value& V = args.at(key);
     if (!V) ExitWithMsg("missing arg values for", key);
     try {std::string tvar;
         tvar = V.asString();
         var = std::stof(tvar);}
     catch (...){ ExitWithMsg("Expected double/float for", key);}
}
void process_arg(std::string& var, argmap& args, std::string key, stag_){
     docopt::value& V = args.at(key);
     if (!V) ExitWithMsg("missing arg values", key);
     try {var = V.asString();} catch (...) { ExitWithMsg("Expected string for", key);}
}
void process_arg(bool& var, argmap& args, std::string key, stag_){
     docopt::value& V = args.at(key);
     if (!V) ExitWithMsg("missing arg values", key);
     try {var = V.asBool();} catch (...) { ExitWithMsg("Expected bool for", key);}
}

#endif
