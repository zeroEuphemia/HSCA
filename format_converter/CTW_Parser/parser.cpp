#include <string>
#include <vector>
#include <cstring>
#include <cstdlib>
#include <unordered_map>
#include <iostream>
#include <fstream>
using namespace std;

enum ParamType { PT_Bool, PT_Enum, PT_Range };

class Param{
public:
    ParamType pt;
    string name;
    virtual ~Param(){}
};

class BoolParam: public Param{
public:
    BoolParam(){
        pt = PT_Bool;
    }
    ~BoolParam(){}
};

class EnumParam: public Param{
public:
    vector<string> values;
    EnumParam(const vector<string>& _v): values(_v){
        pt = PT_Enum;
    }
    ~EnumParam(){}
};

class RangeParam: public Param{
public:
    int l, r;
    int step;
    RangeParam(int _l, int _r, int _step): l(_l), r(_r), step(_step){
        pt = PT_Range;
    }
    ~RangeParam(){}
};

void print_param(const Param* p){
    // cout << "param name: " << p->name << endl;
    // switch (p->pt){
    // case PT_Bool: {
    //     cout << "param type: bool" << endl;
    // }
    //     break;
    
    // case PT_Enum: {
    //     cout << "param type: enum (values: ";
    //     const EnumParam* ep = dynamic_cast<const EnumParam*>(p);
    //     for (const string& v: ep->values){
    //         cout << v << ";";
    //     }
    //     cout << ")" << endl;
    // }
    //     break;
    
    // case PT_Range: {
    //     const RangeParam* rp = dynamic_cast<const RangeParam*>(p);
    //     cout << "param type: range (from " << rp->l << " to " << rp->r << " with step=" << rp->step << ")" << endl;
    // }
    //     break;
    // }
}

enum ImplOperator { IO_Impl, IO_Iff };
enum EqOperator { EO_Eq, EO_Ne };
enum RelOperator { RO_Ge, RO_Gt, RO_Le, RO_Lt };
enum PMOperator { PMO_Plus, PMO_Minus };
enum MMDOperator { MMDO_Mod, MMDO_Mult, MMDO_Div };

enum OrOperator { OO_Unit };
enum AndOperator { AO_Unit };

string tostring_implopr(ImplOperator o){
    switch (o){
    case IO_Impl: return "->";
    case IO_Iff: return "<->";
    }
    return "";
}

string tostring_eqopr(EqOperator o){
    switch (o){
    case EO_Eq: return "==";
    case EO_Ne: return "!=";
    }
    return "";
}

string tostring_relopr(RelOperator o){
    switch (o){
    case RO_Ge: return ">=";
    case RO_Le: return "<=";
    case RO_Gt: return ">";
    case RO_Lt: return "<";
    }
    return "";
}

string tostring_pmopr(PMOperator o){
    switch (o){
    case PMO_Plus: return "+";
    case PMO_Minus: return "-";
    }
    return "";
}

string tostring_mmdopr(MMDOperator o){
    switch (o){
    case MMDO_Div: return "/";
    case MMDO_Mod: return "mod";
    case MMDO_Mult: return "*";
    }
    return "";
}

string tostring_andopr(AndOperator o){
    return "and";
}

string tostring_oropr(OrOperator o){
    return "or";
}

const vector<pair<string, ImplOperator> > dict_implopr{
    {"<=>", IO_Iff}, {"<->", IO_Iff}, {"=>", IO_Impl}, {"->", IO_Impl}, 
};
const vector<pair<string, EqOperator> > dict_eqopr{
    {"==", EO_Eq}, {"=", EO_Eq}, {"!=", EO_Ne}, 
};
const vector<pair<string, RelOperator> > dict_relopr{
    {"<=", RO_Le}, {">=", RO_Ge}, {"<", RO_Lt}, {">", RO_Gt}, 
};
const vector<pair<string, PMOperator> > dict_pmopr{
    {"+", PMO_Plus}, {"-", PMO_Minus}, 
};
const vector<pair<string, MMDOperator> > dict_mmdopr{
    {"%", MMDO_Mod}, {"/", MMDO_Div}, {"*", MMDO_Mult}, 
};
const vector<pair<string, AndOperator> > dict_andopr{
    {"&&", AO_Unit}, {"&", AO_Unit}, {"and", AO_Unit}, {"AND", AO_Unit}, 
};
const vector<pair<string, OrOperator> > dict_oropr{
    {"||", OO_Unit}, {"|", OO_Unit}, {"or", OO_Unit}, {"OR", OO_Unit}, 
};

const vector<string> all_opr_string{
    "<=>", "<->", "=>", "->", "==", "=", "!=", "<=", ">=", "<", ">", 
    "+", "-", "%", "/", "*", "&&", "&", "||", "|", 
    "and", "AND", "or", "OR"
};

template <class Operand, class Operator>
struct BinOp{
    vector<Operand*> opr;
    vector<Operator> opt;

    ~BinOp(){
        for (Operand* p: opr){
            delete p;
        }
    }
};

bool strict_printing;

template <class Operand, class Operator>
void print_binop(const BinOp<Operand, Operator>* p, int spacewidth, string binop_name, string (*tostring_opr)(Operator), void (*print_operand)(const Operand*, int), bool goodcheck){
    int sz = p->opt.size();

    if (goodcheck){
        // if (sz > 1){
        //     cout << "FATAL ERROR! shape check failing at " << binop_name << endl;
        // }
    }
    
    // if (sz > 0 || strict_printing){
    //     for (int i = 0; i < spacewidth; ++i)
    //         cout << ' ';
    //     cout << binop_name << ":" << endl;
    // }
    
    for (int i = 0; i < sz; ++i){
        (*print_operand)(p->opr[i], spacewidth + 2);
        // for (int j = 0; j < spacewidth + 2; ++j)
        //     cout << ' ';
        // cout << (*tostring_opr)(p->opt[i]) << endl;
    }

    (*print_operand)(p->opr[sz], spacewidth + ((sz == 0 && !strict_printing) ? 0: 2));
}

enum SingleType { ST_Not, ST_Impl, ST_Atom };

class Primary{
public:
    SingleType pt; 
    virtual ~Primary(){}
};

typedef BinOp<Primary, MMDOperator> ModMultDiv;
typedef BinOp<ModMultDiv, PMOperator> PlusMinus;
typedef BinOp<PlusMinus, RelOperator> RelExpr;
typedef BinOp<RelExpr, EqOperator> EqExpr;
typedef BinOp<EqExpr, AndOperator> AndExpr;
typedef BinOp<AndExpr, OrOperator> OrExpr;
typedef BinOp<OrExpr, ImplOperator> ImplExpr;


class NotExpr: public Primary{
public:
    Primary* inner;
    NotExpr(Primary* _i): inner(_i){
        pt = ST_Not;
    }
    ~NotExpr(){
        delete inner;
    }
};

class PrimaryImplExpr: public Primary{
public:
    ImplExpr* inner;
    PrimaryImplExpr(ImplExpr* _i): inner(_i){
        pt = ST_Impl;
    }
    ~PrimaryImplExpr(){
        delete inner;
    }
};

enum AtomType { AT_Bool, AT_Value, AT_Int };

class AtomExpr: public Primary{
public:
    AtomType at;
    virtual ~AtomExpr(){}
};

class BoolAtom: public AtomExpr{
public:
    bool val;
    BoolAtom(bool _v): val(_v){
        pt = ST_Atom;
        at = AT_Bool;
    }
    ~BoolAtom(){}
};


class ValueAtom: public AtomExpr{
public:
    string val;
    ValueAtom(string _v): val(_v){
        pt = ST_Atom;
        at = AT_Value;
    }
    ~ValueAtom(){}
};

class IntAtom: public AtomExpr{
public:
    int val;
    IntAtom(int _v): val(_v){
        pt = ST_Atom;
        at = AT_Int;
    }
    ~IntAtom(){}
};

void print_implexpr(const ImplExpr* p, int sw);

void print_atomexpr(const AtomExpr* p, int sw, bool print_space){
    if (print_space){
        // for (int j = 0; j < sw; ++j)
        //     cout << ' ';
    }

    switch (p->at){
    case AT_Bool: {
        // cout << "Boolean: ";
        // const BoolAtom* ba = dynamic_cast<const BoolAtom*>(p);
        // if (ba->val) cout << "TRUE";
        // else cout << "FALSE";
    }
        break;
    
    case AT_Value: {
        // cout << "Value: ";
        // const ValueAtom* va = dynamic_cast<const ValueAtom*>(p);
        // cout << va->val;
    }
        break;
    
    case AT_Int: {
        // cout << "Int: ";
        // const IntAtom* ia = dynamic_cast<const IntAtom*>(p);
        // cout << ia->val;
    }
        break;
    }
    // cout << endl;
}

void print_primary(const Primary* p, int sw){
    switch (p->pt){
    case ST_Not: {
        // for (int j = 0; j < sw; ++j)
        //     cout << ' ';
        // cout << "NotExpression:" << endl;
        const NotExpr* ne = dynamic_cast<const NotExpr*>(p);
        print_primary(ne->inner, sw + 2);
    }
        break;
    
    case ST_Impl: {
        const PrimaryImplExpr* pie = dynamic_cast<const PrimaryImplExpr*>(p);
        print_implexpr(pie->inner, sw);
    }
        break;
    
    case ST_Atom: {
        const AtomExpr* ae = dynamic_cast<const AtomExpr*>(p);
        print_atomexpr(ae, sw, true);
    }
        break;
    }
}

void print_modmultdiv(const ModMultDiv* p, int sw){
    print_binop(p, sw, "ModMultDiv", tostring_mmdopr, print_primary, false);
}

void print_plusminus(const PlusMinus* p, int sw){
    print_binop(p, sw, "PlusMinus", tostring_pmopr, print_modmultdiv, false);
}

void print_relexpr(const RelExpr* p, int sw){
    print_binop(p, sw, "RelationalExpression", tostring_relopr, print_plusminus, true);
}

void print_eqexpr(const EqExpr* p, int sw){
    print_binop(p, sw, "EquationalExpression", tostring_eqopr, print_relexpr, true);
}

void print_andexpr(const AndExpr* p, int sw){
    print_binop(p, sw, "AndExpression", tostring_andopr, print_eqexpr, false);
}

void print_orexpr(const OrExpr* p, int sw){
    print_binop(p, sw, "OrExpression", tostring_oropr, print_andexpr, false);
}

void print_implexpr(const ImplExpr* p, int sw){
    print_binop(p, sw, "ImplicationExpression", tostring_implopr, print_orexpr, true);
}

struct CITModel{
    string model_name;
    vector<Param*> params;
    vector<ImplExpr*> constraints;

    string text;
    int textlen;

    bool ok;
    string error_msg;
    int fail_position;

    void read_first(const char* fname){
        static char tmp[100005];
        FILE* fin = fopen(fname, "r");
        while (fgets(tmp, 100001, fin)){
            int st = 0;
            int len = strlen(tmp);
            while (st < len && isspace(tmp[st]))
                ++st;
            if (tmp[st] == '/' && tmp[st + 1] == '/'){
                continue;
            }
            while (st < len){
                while (st < len && isspace(tmp[st]))
                    ++st;
                int st2 = st;
                while (st2 < len && !isspace(tmp[st2]))
                    ++st2;
                if (st2 > st){
                    
                    if (!text.empty()) text.push_back(' ');
                    while (st < st2){
                        text.push_back(tmp[st]);
                        ++st;
                    }
                }
            }
        }
        fclose(fin);
    }

    void read_spaces(int& st){
        while (st < textlen && isspace(text[st])) ++st;
    }

    void skip_one_empty(int& st){
        if (st < textlen){
            if (text[st] != ' '){
                ok = false;
                error_msg = "no good empty separator";
                fail_position = st;
                return ;
            }
            ++st;
        }
    }

    bool allowed_id_first_char(char c){
        return (c >= 'a' && c <= 'z') || (c >= 'A' && c <= 'Z') || c == '$' || c == '_';
    }

    bool allowed_id_char(char c){
        return (c >= '0' && c <= '9') || allowed_id_first_char(c);
    }

    string read_id(int& st){
        string tmps;
        read_spaces(st);
        if (st >= textlen){
            ok = false;
            error_msg = "empty id";
            fail_position = st;
            return "";
        } else if (!allowed_id_first_char(text[st])){
            ok = false;
            error_msg = "not allowed first char";
            fail_position = st;
            return "";
        }

        do {
            tmps.push_back(text[st]);
            ++st;   
        } while (st < textlen && allowed_id_char(text[st]));
        return tmps;
    }

    void find_next_char_must(char c, int& st){
        read_spaces(st);
        if (st >= textlen || text[st] != c){
            ok = false;
            error_msg = "failed to get char ";
            error_msg.push_back(c);
            fail_position = st;
            return ;
        }
        ++st;
    }

    bool find_next_char(char c, int& st){
        read_spaces(st);
        if (st < textlen && text[st] == c){
            ++st;
            return true;
        }
        return false;
    }

    char read_next_char(int& st){
        read_spaces(st);
        if (st >= textlen){
            ok = false;
            error_msg = "no more char";
            fail_position = st;
            return '!';
        }
        ++st;
        return text[st - 1];
    }

    int read_int(int& st){
        int oldst = st;
        int x = 0;
        while (st < textlen && text[st] >= '0' && text[st] <= '9'){
            x = x * 10 + text[st] - '0';
            ++st;
        }
        if (st == oldst){
            ok = false;
            error_msg = "failed to read int";
            fail_position = st;
            return -1;
        }
        return x;
    }

    bool check_startwith(const string& s, int& st){
        int l = s.length();
        if (st + l > textlen)
            return false;
        for (int i = 0; i < l; ++i)
            if (s[i] != text[st + i])
                return false;
        st += l;
        return true;
    }

    int read_possibly_signed_int(int& st){
        int f = 1;
        if (find_next_char('-', st)){
            f = -1;
        }
        f *= read_int(st);
        if (!ok) return -1;
        return f;
    }

    void find_separator(char c, int& st){
        int oldst = st;
        find_next_char_must(c, st);
        if (ok) return ;

        ok = true;
        error_msg.clear();
        st = oldst;
        skip_one_empty(st);
    }

    void read_range_param(int& st){
        int l = read_possibly_signed_int(st);
        if (!ok) return ;

        char c = read_next_char(st);
        if (!ok) return ;
        if (c != '.'){
            ok = false;
            error_msg = "failed to read the first dot";
            fail_position = st;
            return ;
        }

        if (st >= textlen || text[st] != '.'){
            ok = false;
            error_msg = "failed to read the second dot";
            fail_position = st;
            return ;
        }
        ++st;

        int r = read_possibly_signed_int(st);
        if (!ok) return ;

        find_next_char_must(']', st);
        if (!ok) return ;

        int oldst = st;
        int step = 1;
        string tmps = read_id(st);
        if (!ok){
            st = oldst;
            ok = true;
            error_msg.clear();
        } else {
            if (tmps == "step"){
                step = read_possibly_signed_int(st);
                if (!ok){
                    error_msg = "failed to read step int";
                    fail_position = st;
                    return ;
                }
                if (step < 0){
                    ok = false;
                    error_msg = "step int < 0";
                    fail_position = st;
                    return ;
                }
            } else {
                st = oldst;
            }
        }

        Param* rp = new RangeParam(l, r, step);
        params.push_back(rp);
    }

    void read_boolean_param(char fc, int& st){
        find_separator(',', st);
        if (!ok) return ;

        string next_elem = read_id(st);
        if (!ok) return ;
        if ((fc == 'T' && next_elem != "FALSE") ||
            (fc == 'F' && next_elem != "TRUE") ||
            (fc == 't' && next_elem != "false") ||
            (fc == 'f' && next_elem != "true")){
            ok = false;
            error_msg = "Boolean const not complementary";
            fail_position = st;
            return ;
        }

        find_next_char_must('}', st);
        if (!ok) return ;
    
        Param* bp = new BoolParam();
        params.push_back(bp);
    }

    void read_enum_param(int& st){
        vector<string> values;
        for (; ; ){
            int oldst = st;
            string next_elem = read_id(st);
            if (!ok){
                ok = true;
                error_msg.clear();
                st = oldst;
                find_next_char_must('}', st);
                if (!ok){
                    error_msg = "enum values ending incorrectly";
                    return ;
                }
                break;
            }
            
            values.emplace_back(next_elem);

            oldst = st;
            find_separator(',', st);
            if (!ok){
                ok = true;
                error_msg.clear();
                st = oldst;
                find_next_char_must('}', st);
                if (!ok){
                    error_msg = "enum values ending incorrectly";
                    return ;
                }
                break;
            }
        }

        Param* ep = new EnumParam(values);
        params.push_back(ep);
    }

    void read_model_name(int& st){
        string model_verbatim = read_id(st);
        if (!ok) return;
        if (model_verbatim != "Model"){
            ok = false;
            error_msg = "Model not found";
            fail_position = st;
            return ;
        }

        skip_one_empty(st);
        if (!ok) return ;
        model_name = read_id(st);
    }

    void read_parameters(int& st){
        string parameters_verbatim = read_id(st);
        if (!ok) return ;
        if (parameters_verbatim != "Parameters"){
            ok = false;
            error_msg = "Parameters not found";
            fail_position = st;
            return ;
        }

        find_next_char_must(':', st);
        if (!ok) return;

        for (; ; ){
            int oldst = st;
            string next_id = read_id(st);
            if (!ok){
                if (error_msg[0] == 'e'){
                    ok = true;
                    error_msg.clear();
                    st = oldst;
                    break;
                } else return ;
            }
            if (next_id == "Constraints"){
                st = oldst;
                break ;
            }

            find_next_char_must(':', st);
            if (!ok) return;

            oldst = st;
            char c = read_next_char(st);
            if (!ok) return ;
            if (c == '{'){
                oldst = st;
                string next_elem = read_id(st);
                if (!ok) return ;
                if (next_elem == "TRUE" || next_elem == "FALSE" || next_elem == "true" || next_elem == "false"){
                    read_boolean_param(next_elem[0], st);
                } else {
                    st = oldst;
                    read_enum_param(st);
                }
                if (!ok) return ;
            } else if (c == '['){
                read_range_param(st);
                if (!ok) return ;
            } else if (c == 'B'){
                st = oldst;
                string boolean_verbatim = read_id(st);
                if (!ok) return;
                if (boolean_verbatim != "Boolean"){
                    ok = false;
                    error_msg = "read Boolean failed";
                    fail_position = st;
                    return ;
                }
                
                Param* bp = new BoolParam();
                params.push_back(bp);
            } else {
                ok = false;
                error_msg = "wrong type of param";
                fail_position = st;
                return ;
            }

            params.back()->name = next_id;

            find_separator(';', st);
            if (!ok) return ;
        }

        if (params.empty()){
            ok = false;
            error_msg = "no param found";
            fail_position = st;
            return ;
        }
    }

    template <class Operand, class Operator>
    BinOp<Operand, Operator>* read_binop(const vector<pair<string, Operator> >& dict, Operand* (CITModel::*frec)(int&), int& st){
        BinOp<Operand, Operator>* res = new BinOp<Operand, Operator>();

        for (; ; ){
            Operand* o1 = (this->*frec)(st);
            if (!ok){
                if (o1 != nullptr) delete o1; 
                delete res;
                return nullptr;
            }

            res->opr.push_back(o1);

            bool found = false;
            read_spaces(st);
            for (const pair<string, Operator>& pp: dict){
                int oldst = st;
                if (check_startwith(pp.first, st)){
                    int newst = st;
                    bool noother = true;
                    for (const string& s: all_opr_string){
                        st = oldst;
                        if (s.length() > pp.first.length() && check_startwith(s, st)){
                            noother = false;
                            break;
                        }
                    }
                    if (!noother){
                        st = oldst;
                        continue;
                    }

                    st = newst;
                    res->opt.push_back(pp.second);
                    found = true;
                    break;
                }
            }

            if (!found) break; 
        }
        
        return res;        
    }

    Primary* read_primary(int& st){
        int oldst = st;
        char c = read_next_char(st);
        if (!ok) return nullptr;

        if (c == '('){
            ImplExpr* ie = read_implexpr(st);
            if (!ok) return nullptr;

            find_next_char_must(')', st);
            if (!ok){
                delete ie;
                return nullptr;
            }

            Primary* res = new PrimaryImplExpr(ie);
            return res;
        } else if (c == '!'){
            Primary* p = read_primary(st);
            if (!ok) return nullptr;

            Primary* res = new NotExpr(p);
            return res;
        } else {
            st = oldst;
            string next_elem = read_id(st);
            if (!ok){
                if (error_msg[0] == 'e') return nullptr;

                ok = true;
                error_msg.clear();
                st = oldst;
                int x = read_possibly_signed_int(st);
                if (!ok){
                    error_msg = "cannot find integer atom";
                    return nullptr;
                }

                Primary* res = new IntAtom(x);
                return res;
            }
            if (next_elem == "not"){
                skip_one_empty(st);
                if (!ok) return nullptr;

                Primary* p = read_primary(st);
                if (!ok) return nullptr;

                Primary* res = new NotExpr(p);
                return res;
            }
            if (next_elem == "true" || next_elem == "TRUE"){
                Primary* res = new BoolAtom(true);
                return res;
            }
            if (next_elem == "false" || next_elem == "FALSE"){
                Primary* res = new BoolAtom(false);
                return res;
            }
            Primary* res = new ValueAtom(next_elem);
            return res;
        }

        return nullptr;
    }

    ModMultDiv* read_modmultdiv(int& st){
        return read_binop(dict_mmdopr, &CITModel::read_primary, st);
    }

    PlusMinus* read_plusminus(int& st){
        return read_binop(dict_pmopr, &CITModel::read_modmultdiv, st);
    }

    RelExpr* read_relexpr(int& st){
        return read_binop(dict_relopr, &CITModel::read_plusminus, st);
    }

    EqExpr* read_eqexpr(int& st){
        return read_binop(dict_eqopr, &CITModel::read_relexpr, st);
    }

    AndExpr* read_andexpr(int& st){
        return read_binop(dict_andopr, &CITModel::read_eqexpr, st);
    }

    OrExpr* read_orexpr(int& st){
        return read_binop(dict_oropr, &CITModel::read_andexpr, st);
    }

    ImplExpr* read_implexpr(int& st){
        return read_binop(dict_implopr, &CITModel::read_orexpr, st);
    }

    void read_constraints(int& st){
        string constraints_verbatim = read_id(st);
        if (!ok){
            if (error_msg[0] == 'e'){
                ok = true;
                error_msg.clear();
            }
            return ;
        }
        if (constraints_verbatim != "Constraints"){
            ok = false;
            error_msg = "Constraints not found";
            fail_position = st;
            return ;
        }

        find_next_char_must(':', st);
        if (!ok) return;

        for (; ; ){
            find_next_char_must('#', st);
            if (!ok){
                if (st >= textlen){
                    ok = true;
                    error_msg.clear();
                    break;
                }
                return ;
            }

            ImplExpr* res = read_implexpr(st);
            if (!ok) return ;

            constraints.push_back(res);

            find_next_char_must('#', st);
            if (!ok) return;
        }
    }

    void read(const char* fname){
        read_first(fname);
        textlen = text.length();

        int st = 0;
        ok = true;
        error_msg = "";

        read_model_name(st);
        if (!ok) return ;

        skip_one_empty(st);
        if (!ok) return ;
        read_parameters(st);
        if (!ok) return ;

        read_constraints(st);
    }

    void print_model() const {
        // cout << "model name: " << model_name << endl;
        // cout << endl;

        for (const Param* p: params){
            print_param(p);
            // cout << endl;
        }

        for (const ImplExpr* ie: constraints){
            print_implexpr(ie, 0);
            // cout << endl;
        }
    }

    unordered_map<string, const Param*> dump_param_info(){
        unordered_map<string, const Param*> res;
        for (const Param* p: params){
            res[p->name] = p;
        }
        return res;
    }

    vector<string> dump_param_names(){
        vector<string> res;
        for (const Param* p: params){
            res.push_back(p->name);
        }
        return res;
    }

    ~CITModel(){
        for (Param* p: params) delete p;
        for (ImplExpr* p: constraints) delete p;
    }
};


class SimpleNode {
public:
    bool notmark;
    bool isleaf;
    virtual ~SimpleNode(){}
};

class SimpleBranchNode: public SimpleNode {
public:
    vector<string> opt;
    vector<SimpleNode*> opr;
    SimpleBranchNode(){
        notmark = false;
        isleaf = false;
    }
    ~SimpleBranchNode(){
        for (SimpleNode* p: opr){
            delete p;
        }
    }
};

class SimpleLeafNode: public SimpleNode {
public:
    AtomExpr* ae;
    SimpleLeafNode(){
        notmark = false;
        isleaf = true;
        ae = nullptr;
    }
    ~SimpleLeafNode(){ delete ae; }
};

AtomExpr* clone_atom_expr(const AtomExpr* p){
    AtomExpr* res;

    switch (p->at){
    case AT_Bool: {
        const BoolAtom* ba = dynamic_cast<const BoolAtom*>(p);
        res = new BoolAtom(ba->val);
    }
        break;
    
    case AT_Value: {
        const ValueAtom* va = dynamic_cast<const ValueAtom*>(p);
        res = new ValueAtom(va->val);
    }
        break;
    
    case AT_Int: {
        const IntAtom* ia = dynamic_cast<const IntAtom*>(p);
        res = new IntAtom(ia->val);
    }
        break;
    }

    return res;
}

SimpleNode* clone_simple_node(const SimpleNode* p){
    if (p->isleaf){
        const SimpleLeafNode* sln = dynamic_cast<const SimpleLeafNode*>(p);
        SimpleLeafNode* res = new SimpleLeafNode();
        res->notmark = sln->notmark;
        res->ae = clone_atom_expr(sln->ae);
        return res;
    }

    const SimpleBranchNode* sbn = dynamic_cast<const SimpleBranchNode*>(p);
    SimpleBranchNode* res = new SimpleBranchNode();
    res->notmark = sbn->notmark;
    res->opt = sbn->opt;
    for (const SimpleNode* sn: sbn->opr){
        res->opr.push_back(clone_simple_node(sn));
    }
    return res;
}

template <class Operand, class Operator>
SimpleNode* compile_binop(const BinOp<Operand, Operator>* p, string (*tostring_opr)(Operator), SimpleNode* (*compile_operand)(const Operand*)){
    int sz = p->opt.size();
    if (sz == 0) return compile_operand(p->opr[0]); 

    SimpleBranchNode* sbn = new SimpleBranchNode();
    
    for (int i = 0; i < sz; ++i){
        sbn->opt.push_back((*tostring_opr)(p->opt[i]));
        sbn->opr.push_back(compile_operand(p->opr[i]));
    }

    sbn->opr.push_back(compile_operand(p->opr[sz]));
    return sbn;
}

SimpleNode* compile_implexpr(const ImplExpr* p);

SimpleNode* compile_atomexpr(const AtomExpr* p){
    SimpleLeafNode* res = new SimpleLeafNode();

    switch (p->at){
    case AT_Bool: {
        const BoolAtom* ba = dynamic_cast<const BoolAtom*>(p);
        BoolAtom* ba2 = new BoolAtom(ba->val);
        res->ae = ba2;
    }
        break;
    
    case AT_Value: {
        const ValueAtom* va = dynamic_cast<const ValueAtom*>(p);
        ValueAtom* va2 = new ValueAtom(va->val);
        res->ae = va2;
    }
        break;
    
    case AT_Int: {
        const IntAtom* ia = dynamic_cast<const IntAtom*>(p);
        IntAtom* ia2 = new IntAtom(ia->val);
        res->ae = ia2;
    }
        break;
    }

    return res;
}

SimpleNode* compile_primary(const Primary* p){
    SimpleNode* res = nullptr;

    switch (p->pt){
    case ST_Not: {
        const NotExpr* ne = dynamic_cast<const NotExpr*>(p);
        res = compile_primary(ne->inner);
        res->notmark = !res->notmark;
    }
        break;
    
    case ST_Impl: {
        const PrimaryImplExpr* pie = dynamic_cast<const PrimaryImplExpr*>(p);
        res = compile_implexpr(pie->inner);
    }
        break;
    
    case ST_Atom: {
        const AtomExpr* ae = dynamic_cast<const AtomExpr*>(p);
        res = compile_atomexpr(ae);
    }
        break;
    }

    return res;
}

SimpleNode* compile_modmultdiv(const ModMultDiv* p){
    return compile_binop(p, tostring_mmdopr, compile_primary);
}

SimpleNode* compile_plusminus(const PlusMinus* p){
    return compile_binop(p, tostring_pmopr, compile_modmultdiv);
}

SimpleNode* compile_relexpr(const RelExpr* p){
    return compile_binop(p, tostring_relopr, compile_plusminus);
}

SimpleNode* compile_eqexpr(const EqExpr* p){
    return compile_binop(p, tostring_eqopr, compile_relexpr);
}

SimpleNode* compile_andexpr(const AndExpr* p){
    return compile_binop(p, tostring_andopr, compile_eqexpr);
}

SimpleNode* compile_orexpr(const OrExpr* p){
    return compile_binop(p, tostring_oropr, compile_andexpr);
}

SimpleNode* compile_implexpr(const ImplExpr* p){
    return compile_binop(p, tostring_implopr, compile_orexpr);
}

void print_simple_node(const SimpleNode* p, int sw){
    // for (int i = 0; i < sw; ++i)
    //     cout << ' ';
    // if (p->notmark) cout << "not ";

    if (p->isleaf){
        const SimpleLeafNode* sln = dynamic_cast<const SimpleLeafNode*>(p);
        print_atomexpr(sln->ae, sw, false);
        return ;
    }

    const SimpleBranchNode* sbn = dynamic_cast<const SimpleBranchNode*>(p);
    int sz = (int) sbn->opt.size();
    // cout << "(" << endl;

    // if (sz == 0){ cout << "FATAL ERROR! arity wrong!" << endl; }
    
    for (int i = 0; i < sz; ++i){
        print_simple_node(sbn->opr[i], sw + 2);
        // for (int j = 0; j < sw + 2; ++j)
        //     cout << ' ';
        // cout << sbn->opt[i] << endl;
    }

    print_simple_node(sbn->opr[sz], sw + 2);
    // for (int i = 0; i < sw; ++i)
    //     cout << ' ';
    // cout << ")" << endl;
}


void print_simple_node_simple(const SimpleNode* p, int sw){
    // for (int i = 0; i < sw; ++i)
    //     cout << ' ';
    // if (p->notmark) cout << "not ";

    if (p->isleaf){
        const SimpleLeafNode* sln = dynamic_cast<const SimpleLeafNode*>(p);
        print_atomexpr(sln->ae, sw, false);
        return ;
    }

    const SimpleBranchNode* sbn = dynamic_cast<const SimpleBranchNode*>(p);
    int sz = (int) sbn->opr.size();
    // cout << sbn->opt[0] << " (" << endl;

    // if (sz == 1){ cout << "FATAL ERROR! arity wrong!" << endl; }
    
    for (int i = 0; i < sz; ++i){
        print_simple_node_simple(sbn->opr[i], sw + 2);
    }

    // for (int i = 0; i < sw; ++i)
    //     cout << ' ';
    // cout << ")" << endl;
}

string opt_negate(const string& s){
    if (s == ">="){
        return "<";
    } else if (s == "<="){
        return ">";
    } else if (s == ">"){
        return "<=";
    } else if (s == "<"){
        return ">=";
    }
    string res = s;
    res[0] = '!' + '=' - s[0];
    return res;
}


SimpleNode* normalize(SimpleNode* p, bool relop_passed, const unordered_map<string, const Param*>& pmp){
    if (p->isleaf){
        SimpleLeafNode* sln = dynamic_cast<SimpleLeafNode*>(p);
        if (!relop_passed){
            if (sln->ae->at != AT_Value){
                // cout << "FATAL ERROR! single atom is not param type!" << endl;
                return p;
            }

            const ValueAtom* va = dynamic_cast<const ValueAtom*>(sln->ae);
            if (!pmp.count(va->val)){
                // cout << "FATAL ERROR! single atom is not a param!" << endl;
                return p;
            }

            const Param* param = pmp.find(va->val)->second;
            if (param->pt != PT_Bool){
                // cout << "FATAL ERROR! single atom is not a bool param!" << endl;
                return p;
            }


            SimpleLeafNode* slnb = new SimpleLeafNode();
            slnb->ae = new BoolAtom(!sln->notmark);
            sln->notmark = false;

            SimpleBranchNode* sbn = new SimpleBranchNode();
            sbn->opt.push_back("==");
            sbn->opr.push_back(sln);
            sbn->opr.push_back(slnb);

            p = sbn;
        }

        return p;
    }

    SimpleBranchNode* sbn = dynamic_cast<SimpleBranchNode*>(p);
    int sz = (int) sbn->opt.size();

    bool thisisrel = (sz == 1 && (sbn->opt[0] == "==" ||
                                sbn->opt[0] == "!=" ||
                                sbn->opt[0] == "<=" ||
                                sbn->opt[0] == ">=" ||
                                sbn->opt[0] == "<" ||
                                sbn->opt[0] == ">"));
    
    for (int i = 0; i <= sz; ++i){
        sbn->opr[i] = normalize(sbn->opr[i], thisisrel, pmp);
    }

    if (thisisrel){
        if (!sbn->opr[0]->isleaf || !sbn->opr[1]->isleaf){
            // cout << "FATAL ERROR! one side of rel is not leaf!" << endl;
            return sbn;
        }

        SimpleLeafNode* slnl = dynamic_cast<SimpleLeafNode*>(sbn->opr[0]);
        if (slnl->ae->at != AT_Value){
            // cout << "FATAL ERROR! left side is not param type!" << endl;
            return sbn;
        }
        
        const ValueAtom* va = dynamic_cast<const ValueAtom*>(slnl->ae);
        if (!pmp.count(va->val)){
            // cout << "FATAL ERROR! left side is not a param!" << endl;
            return sbn;
        }

        const Param* param = pmp.find(va->val)->second;
        SimpleLeafNode* slnr = dynamic_cast<SimpleLeafNode*>(sbn->opr[1]);

        if (sbn->notmark){
            sbn->opt[0] = opt_negate(sbn->opt[0]);
            sbn->notmark = false;
        }

        switch (param->pt){
        case PT_Bool: {
            if (slnr->ae->at != AT_Bool){
                // cout << "FATAL ERROR! left side is bool, right side is not" << endl;
                return sbn;
            }

            if (sbn->opt[0] != "==" && sbn->opt[0] != "!="){
                // cout << "FATAL ERROR! operator does not fit bool" << endl;
                return sbn;
            }

            bool change = !((slnl->notmark && slnr->notmark) || (!slnl->notmark && !slnr->notmark));
            if (change){
                sbn->opt[0][0] = '!' + '=' - sbn->opt[0][0];
            }
            slnl->notmark = slnr->notmark = false;
        }
            break;
        
        case PT_Enum: {
            if (slnr->ae->at != AT_Value){
                // cout << "FATAL ERROR! left side is enum, right side is not value" << endl;
                return sbn;
            }
            const EnumParam* ep = dynamic_cast<const EnumParam*>(param);
            const ValueAtom* var = dynamic_cast<const ValueAtom*>(slnr->ae);
            bool ok = false;
            for (const string& s: ep->values){
                if (s == var->val){
                    ok = true;
                    break;
                }
            }

            if (!ok){
                // cout << "FATAL ERROR! left side is enum, cannot take right side value" << endl;
                return sbn;
            }
        }
            break;
        
        case PT_Range: {
            if (slnr->ae->at != AT_Int){
                // cout << "FATAL ERROR! left side is range, right side is not int" << endl;
                return sbn;
            }
        }
            break;
        }
    }

    return sbn;
}


SimpleNode* compile_iff_impl(SimpleNode* p){
    if (p->isleaf) return p;

    SimpleBranchNode* sbn = dynamic_cast<SimpleBranchNode*>(p);
    int sz = (int) sbn->opt.size();

    bool thisisrel = (sz == 1 && (sbn->opt[0] == "==" ||
                                sbn->opt[0] == "!=" ||
                                sbn->opt[0] == "<=" ||
                                sbn->opt[0] == ">=" ||
                                sbn->opt[0] == "<" ||
                                sbn->opt[0] == ">"));
    
    if (thisisrel){
        if (sbn->notmark){
            sbn->opt[0] = opt_negate(sbn->opt[0]);
            sbn->notmark = false;
        }

        return sbn;
    }

    if (sbn->opt[0] == "->"){
        sbn->opt[0] = "or";
        sbn->opr[0]->notmark = !sbn->opr[0]->notmark;
    } else if (sbn->opt[0] == "<->"){
        SimpleNode* slnl2 = clone_simple_node(sbn->opr[0]);
        SimpleNode* slnr2 = clone_simple_node(sbn->opr[1]);
        sbn->opr[0]->notmark = !sbn->opr[0]->notmark;
        sbn->opr[1]->notmark = !sbn->opr[1]->notmark;
        SimpleBranchNode* sbnnewl = new SimpleBranchNode();
        SimpleBranchNode* sbnnewr = new SimpleBranchNode();
        sbnnewl->opt.push_back("or");
        sbnnewl->opr.push_back(sbn->opr[0]);
        sbnnewl->opr.push_back(slnr2);
        sbnnewr->opt.push_back("or");
        sbnnewr->opr.push_back(sbn->opr[1]);
        sbnnewr->opr.push_back(slnl2);
        sbn->opt[0] = "and";
        sbn->opr[0] = sbnnewl;
        sbn->opr[1] = sbnnewr;
    }

    bool exists_not_or = false, exists_not_and = false;
    for (const string& s: sbn->opt){
        if (s != "or") exists_not_or = true;
        if (s != "and") exists_not_and = true;
    }
    if (exists_not_or && exists_not_and){
        // cout << "FATAL ERROR! unexpected connective" << endl;
        return sbn;
    }

    if (sbn->notmark){
        for (int i = 0; i <= sz; ++i){
            sbn->opr[i]->notmark = !sbn->opr[i]->notmark;
        }
        if (!exists_not_or){
            for (int i = 0; i < sz; ++i){
                sbn->opt[i] = "and";
            }
        } else {
            for (int i = 0; i < sz; ++i){
                sbn->opt[i] = "or";
            }
        }
        sbn->notmark = false;
    }
    
    for (int i = 0; i <= sz; ++i){
        sbn->opr[i] = compile_iff_impl(sbn->opr[i]);
    }

    return sbn;
}


SimpleNode* atomicize(SimpleNode* p, const unordered_map<string, const Param*>& pmp){
    if (p->isleaf) return p;

    SimpleBranchNode* sbn = dynamic_cast<SimpleBranchNode*>(p);
    int sz = (int) sbn->opt.size();

    bool thisisrel = (sz == 1 && (sbn->opt[0] == "==" ||
                                sbn->opt[0] == "!=" ||
                                sbn->opt[0] == "<=" ||
                                sbn->opt[0] == ">=" ||
                                sbn->opt[0] == "<" ||
                                sbn->opt[0] == ">"));
    
    if (thisisrel){
        if (sbn->notmark){
            sbn->opt[0] = opt_negate(sbn->opt[0]);
            sbn->notmark = false;
        }

        SimpleLeafNode* slnl = dynamic_cast<SimpleLeafNode*>(sbn->opr[0]);
        SimpleLeafNode* slnr = dynamic_cast<SimpleLeafNode*>(sbn->opr[1]);
        const ValueAtom* va = dynamic_cast<const ValueAtom*>(slnl->ae);
        const Param* param = pmp.find(va->val)->second;

        if (slnl->notmark || slnr->notmark){
            // cout << "FATAL ERROR! should have no not on opr" << endl;
            return sbn;
        }

        SimpleNode* res;

        switch (param->pt){
        case PT_Bool: {
            const BoolAtom* ba = dynamic_cast<const BoolAtom*>(slnr->ae);
            bool posi = true;
            if (!ba->val) posi = !posi;
            if (sbn->opt[0][0] == '!') posi = !posi;
            
            SimpleLeafNode* slnnew = new SimpleLeafNode();
            slnnew->ae = new ValueAtom(va->val);
            slnnew->notmark = !posi;
            res = slnnew;
        }
            break;
        
        case PT_Enum: {
            const ValueAtom* var = dynamic_cast<const ValueAtom*>(slnr->ae);
            SimpleLeafNode* slnnew = new SimpleLeafNode();
            slnnew->ae = new ValueAtom(va->val + "@" + var->val);
            slnnew->notmark = (sbn->opt[0][0] == '!');
            res = slnnew;
        }
            break;
        
        case PT_Range: {
            const IntAtom* ia = dynamic_cast<const IntAtom*>(slnr->ae);
            const RangeParam* rp = dynamic_cast<const RangeParam*>(param);
            int rl = 1, rr = 0;
            bool posi = true;
            if (sbn->opt[0] == "!="){
                posi = false;
                if (ia->val >= rp->l && ia->val <= rp->r){
                    rl = rr = ia->val;
                }
            } else if (sbn->opt[0] == "=="){
                if (ia->val >= rp->l && ia->val <= rp->r){
                    rl = rr = ia->val;
                }
            } else if (sbn->opt[0][0] == '<'){
                rl = rp->l, rr = ia->val - 1;
                if (sbn->opt[0].length() > 1){
                    ++rr;
                }
                rr = min(rr, rp->r);
            } else if (sbn->opt[0][0] == '>'){
                rl = ia->val + 1, rr = rp->r;
                if (sbn->opt[0].length() > 1){
                    --rl;
                }
                rl = max(rl, rp->l);
            }

            if (rr < rl){
                SimpleBranchNode* sbnnew = new SimpleBranchNode();
                if (posi){
                    sbnnew->opt.push_back("and");
                } else {
                    sbnnew->opt.push_back("or");
                }
                SimpleLeafNode* slnnew = new SimpleLeafNode();
                slnnew->ae = new ValueAtom("@top");
                SimpleLeafNode* slnnew2 = new SimpleLeafNode();
                slnnew2->ae = new ValueAtom("@top");
                slnnew2->notmark = true;

                sbnnew->opr.push_back(slnnew);
                sbnnew->opr.push_back(slnnew2);
                res = sbnnew;
            } else if (rl == rr){
                SimpleLeafNode* slnnew = new SimpleLeafNode();
                slnnew->ae = new ValueAtom(va->val + "@" + to_string(rl));
                slnnew->notmark = !posi;
                res = slnnew;
            } else {
                SimpleBranchNode* sbnnew = new SimpleBranchNode();
                if (posi){
                    for (int i = rl; i < rr; ++i){
                        sbnnew->opt.push_back("or");
                    }
                } else {
                    for (int i = rl; i < rr; ++i){
                        sbnnew->opt.push_back("and");
                    }
                }
                for (int i = rl; i <= rr; ++i){
                    SimpleLeafNode* slnnew = new SimpleLeafNode();
                    slnnew->ae = new ValueAtom(va->val + "@" + to_string(i));
                    slnnew->notmark = !posi;
                    sbnnew->opr.push_back(slnnew);
                }
                res = sbnnew;
            }
        }
            break;
        }

        delete p;
        return res;
    }

    for (int i = 0; i <= sz; ++i){
        sbn->opr[i] = atomicize(sbn->opr[i], pmp);
    }

    return sbn;
}


SimpleNode* flatten(SimpleNode* p){
    if (p->isleaf) return p;

    SimpleBranchNode* sbn = dynamic_cast<SimpleBranchNode*>(p);
    int sz = (int) sbn->opt.size();
    string thisopt = sbn->opt[0];

    vector<SimpleNode*> newvec;
    for (int i = 0; i <= sz; ++i){
        SimpleNode* tmp = sbn->opr[i];
        tmp = flatten(tmp);
        if (!tmp->isleaf){
            SimpleBranchNode* childsbn = dynamic_cast<SimpleBranchNode*>(tmp);
            int childsz = (int) childsbn->opt.size();
            if (childsz == 0){
                // cout << "FATAL ERROR! 0 arity branch node" << endl;
                newvec.push_back(tmp);
                continue ;
            }
            
            if (childsbn->opt[0] != thisopt){
                newvec.push_back(tmp);
                continue;
            }

            for (int j = 0; j <= childsz; ++j){
                newvec.push_back(childsbn->opr[j]);
            }
            childsbn->opr.clear();
            delete childsbn;
        } else {
            newvec.push_back(tmp);
        }
    }

    sbn->opr.clear();
    sbn->opr = newvec;
    sz = (int) newvec.size();
    sbn->opt.clear();
    for (int i = 0; i < sz - 1; ++i){
        sbn->opt.push_back(thisopt);
    }

    return sbn;
}


string get_opt(const SimpleNode* p){
    if (p->isleaf) return "";
    const SimpleBranchNode* sbn = dynamic_cast<const SimpleBranchNode*>(p);
    return sbn->opt[0];
}














    


    













vector<SimpleNode*> encoded;

SimpleLeafNode* tseitin_encode(SimpleNode* p, int& varcount){
    if (p->isleaf){
        SimpleLeafNode* sln = dynamic_cast<SimpleLeafNode*>(p);   
        return sln;
    }

    SimpleBranchNode* sbn = dynamic_cast<SimpleBranchNode*>(p);
    int sz = (int) sbn->opr.size();
    string thisopt = sbn->opt[0];

    vector<SimpleLeafNode*> newvec;
    for (int i = 0; i < sz; ++i){
        newvec.push_back(tseitin_encode(sbn->opr[i], varcount));
    }

    sbn->opt.clear();
    sbn->opr.clear();
    delete sbn;


    ++varcount;
    int thisvar = varcount;
    SimpleLeafNode* newsln = new SimpleLeafNode();
    newsln->ae = new ValueAtom("@tmp" + to_string(thisvar));

    SimpleLeafNode* newsln_copy = new SimpleLeafNode();
    newsln_copy->ae = new ValueAtom("@tmp" + to_string(thisvar));

    if (thisopt == "and"){
        for (SimpleLeafNode* subsln: newvec){
            SimpleBranchNode* newsbn = new SimpleBranchNode();
            newsbn->opt.push_back("or");
            SimpleNode* tmpsln1 = clone_simple_node(newsln);
            SimpleNode* tmpsln2 = clone_simple_node(subsln);
            tmpsln1->notmark = !tmpsln1->notmark;
            newsbn->opr.push_back(tmpsln1);
            newsbn->opr.push_back(tmpsln2);
            encoded.push_back(newsbn);
        }

        SimpleBranchNode* newsbn = new SimpleBranchNode();
        newsbn->opt.push_back("or");
        for (SimpleLeafNode* subsln: newvec){
            subsln->notmark = !subsln->notmark;
            newsbn->opr.push_back(subsln);
        }
        newsbn->opr.push_back(newsln);
        encoded.push_back(newsbn);
    } else {
        for (SimpleLeafNode* subsln: newvec){
            SimpleBranchNode* newsbn = new SimpleBranchNode();
            newsbn->opt.push_back("or");
            SimpleNode* tmpsln1 = clone_simple_node(newsln);
            SimpleNode* tmpsln2 = clone_simple_node(subsln);
            tmpsln2->notmark = !tmpsln2->notmark;
            newsbn->opr.push_back(tmpsln1);
            newsbn->opr.push_back(tmpsln2);
            encoded.push_back(newsbn);
        }

        SimpleBranchNode* newsbn = new SimpleBranchNode();
        newsbn->opt.push_back("or");
        for (SimpleLeafNode* subsln: newvec){
            newsbn->opr.push_back(subsln);
        }
        newsln->notmark = !newsln->notmark;
        newsbn->opr.push_back(newsln);
        encoded.push_back(newsbn);
    }

    return newsln_copy;
}

void process_or_subtree(SimpleNode* p, int& varcount){
    SimpleBranchNode* sbn = dynamic_cast<SimpleBranchNode*>(p);
    int sz = (int) sbn->opr.size();

    for (int i = 0; i < sz; ++i){
        sbn->opr[i] = tseitin_encode(sbn->opr[i], varcount);
    }

    sbn->opt.clear();
    sbn->opt.push_back("or");
    encoded.push_back(sbn);
}



int main(int argc, const char* argv[]){
    if (argc != 4 && argc != 5) return 0;

    strict_printing = (argc == 5);

    CITModel md;
    md.read(argv[1]);

    if (!md.ok){
        // cout << "ERROR" << endl;
        // cout << md.error_msg << endl;
        for (char& c: md.text){
            if (c == ' ')
                c = '~';
        }
        // cout << md.text.substr(0, md.fail_position) << endl;
    } else {
        // cout << "OK" << endl;
        md.print_model();

        const unordered_map<string, const Param*> pmp = md.dump_param_info();

        FILE *fout_param = fopen(argv[3], "w");
        vector<string> param_names = md.dump_param_names();
        for (auto& paname: param_names){
            fprintf(fout_param, "%s\n", paname.c_str());
        }
        fclose(fout_param);

        vector<vector<int> > cnf;


        unordered_map<string, int> var2int;
        int varcount = 0;
        for (auto& kvp: pmp){
            const Param* p = kvp.second;
            switch (p->pt){
            case PT_Bool: {
                ++varcount;
                var2int[p->name] = varcount; 
            }
                break;
            
            case PT_Enum: {
                const EnumParam* ep = dynamic_cast<const EnumParam*>(p);
                const vector<string>& tmpvals = ep->values;
                int oldvarcnt = varcount;
                for (const string& v: tmpvals){
                    ++varcount;
                    var2int[p->name + "@" + v] = varcount;
                }

                vector<int> clause;
                for (int i = oldvarcnt + 1; i <= varcount; ++i){
                    clause.push_back(i);
                }
                cnf.push_back(clause);
                clause.clear();

                for (int i = oldvarcnt + 1; i <= varcount; ++i){
                    for (int j = i + 1; j <= varcount; ++j){
                        clause.push_back(-i);
                        clause.push_back(-j);
                        cnf.push_back(clause);
                        clause.clear();
                    }
                }
            }
                break;
            
            case PT_Range: {
                const RangeParam* rp = dynamic_cast<const RangeParam*>(p);
                int oldvarcnt = varcount;
                for (int i = rp->l; i <= rp->r; ++i){
                    ++varcount;
                    var2int[p->name + "@" + to_string(i)] = varcount;
                }

                vector<int> clause;
                for (int i = oldvarcnt + 1; i <= varcount; ++i){
                    clause.push_back(i);
                }
                cnf.push_back(clause);
                clause.clear();

                for (int i = oldvarcnt + 1; i <= varcount; ++i){
                    for (int j = i + 1; j <= varcount; ++j){
                        clause.push_back(-i);
                        clause.push_back(-j);
                        cnf.push_back(clause);
                        clause.clear();
                    }
                }
            }
                break;
            }
        }
        ++varcount;
        var2int["@top"] = varcount;

        int old_varcount = varcount;

        for (const ImplExpr* ie: md.constraints){
            SimpleNode* sn = compile_implexpr(ie);
            
            sn = normalize(sn, false, pmp);
            sn = compile_iff_impl(sn);
            sn = atomicize(sn, pmp);
            sn = flatten(sn);

            // cout << "\033[;32m prenode \033[0m" << endl;
            print_simple_node(sn, 0);

            string o = get_opt(sn);

            if (o == "or"){
                process_or_subtree(sn, varcount);
            } else if (o == "and"){
                SimpleBranchNode* sbn = dynamic_cast<SimpleBranchNode*>(sn);
                int sz = (int) sbn->opr.size();
                for (int i = 0; i < sz; ++i){
                    SimpleNode* tmp = sbn->opr[i];
                    if (tmp->isleaf){
                        encoded.push_back(tmp);
                    } else {
                        process_or_subtree(tmp, varcount);
                    }
                }
                
                sbn->opr.clear();
                sbn->opt.clear();
                delete sbn;
            } else {
                encoded.push_back(sn);
            }
        }

        for (int vv = old_varcount + 1; vv <= varcount; ++vv){
            var2int["@tmp" + to_string(vv)] = vv;
        }

        cnf.push_back(vector<int>{varcount, -varcount});

        int numcl_original = (int) encoded.size();
        int numcl_enum = (int) cnf.size();

        ofstream cnfout;
        cnfout.open(argv[2]);
        for (auto& kvp: var2int){
            cnfout << "c " << kvp.first << " --> " << kvp.second << endl;
        }
        cnfout << "p cnf " << varcount << " " << (numcl_original + numcl_enum) << endl;
        
        for (int cid = 0; cid < numcl_original; ++cid){
            SimpleNode* sn = encoded[cid];

            if (sn->isleaf){
                SimpleLeafNode* sln = dynamic_cast<SimpleLeafNode*>(sn);
                const ValueAtom* va = dynamic_cast<const ValueAtom*>(sln->ae);
                int l = (sln->notmark ? (-1): 1) * var2int[va->val];
                cnfout << l << " ";
            } else {
                SimpleBranchNode* sbn = dynamic_cast<SimpleBranchNode*>(sn);
                for (SimpleNode* sn2: sbn->opr){
                    SimpleLeafNode* sln = dynamic_cast<SimpleLeafNode*>(sn2);
                    const ValueAtom* va = dynamic_cast<const ValueAtom*>(sln->ae);
                    int l = (sln->notmark ? (-1): 1) * var2int[va->val];
                    cnfout << l << " ";
                }
            }

            cnfout << "0" << endl;

            delete sn;
        }

        for (auto& c: cnf){
            for (auto l: c){
                cnfout << l << " ";
            }
            cnfout << "0" << endl;
        }

        cnfout.close();
    }

    return 0;
}
