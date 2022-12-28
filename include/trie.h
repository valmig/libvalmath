#ifndef TRIE_H_INCLUDED
#define TRIE_H_INCLUDED

#include <string>
#include <initializer_list>

namespace val
{

namespace trie
{
    enum {SIZE_LOWER_SIMPLE = 26, SIZE_CASE_SENSITIVE_SIMPLE = 58};
}

template <class T = std::string, int ALPHA_SIZE = trie::SIZE_LOWER_SIMPLE>
class trie_type
{
private:
    struct trie_node
    {
        trie_node* children[ALPHA_SIZE];
        bool end_of_word = false;
        trie_node() { for (int i = 0; i < ALPHA_SIZE; ++i) children[i] = nullptr; }
    };
    //
    trie_node *head = nullptr;
    int beg = int('a'), nwords = 0;                  // beg = int('a') if valid char starts with 'a' (lower case). For case-sensitive usage,
    //                                                  construct with ALPA_SIZE = 58 and beg = int('A')
    void del(trie_node *node);
    static int valid_word(const T &s, int beg);
    trie_node* Search(const T &prefix) const; // Search node with prefix
    template <template <typename> class C> void appendmatchingword(trie_node *current,C<T> &G,const T &prefix, const T &suffix = T()) const;
    static void copy(trie_node* source, trie_node* &dest);
    //
public:
    trie_type() = default;
    explicit trie_type(int b) : beg(b) {}
    template <template <typename> class C> trie_type(const C<T> &G, int b = int('a')) : beg(b) {insert(G);}
    trie_type(const std::initializer_list<T> &G, int b = int('a')) {for (const auto &v : G) insert(v);}
    ~trie_type() {del(head);}
    trie_type(const trie_type<T,ALPHA_SIZE> &W) : head(nullptr), beg(W.beg), nwords(W.nwords) {copy(W.head,head);}
    trie_type(trie_type<T,ALPHA_SIZE> && Trie) {head = Trie.head; nwords = Trie.nwords; beg = Trie.beg; Trie.head = nullptr; Trie.nwords = 0;}
    //
    const trie_type<T,ALPHA_SIZE>& operator =(const trie_type<T,ALPHA_SIZE>& W);
    const trie_type<T,ALPHA_SIZE>& operator =(trie_type<T,ALPHA_SIZE>&& W);
    //
    void del() {del(head); head = nullptr; nwords = 0;}
    void insert(const T &s);
    template <template <typename> class C> void insert(const C<T> &G) {for (const auto& s : G) insert(s);}
    int numberofelements() const {return nwords;}
    bool containsword(const T& word) const;
    bool containsprefix(const T& prefix) const {return (Search(prefix) != nullptr);}
    template <template <typename> class C> C<T> getmatchingprefix(const T &prefix) const;
};


/*
template <class T, int ALPHA_SIZE>
int trie_type<T,ALPHA_SIZE>::beg = int('a');
*/

template <class T, int ALPHA_SIZE>
int trie_type<T,ALPHA_SIZE>::valid_word(const T &s, int beg)
{
    int l = s.length(), index;
    for (int i = 0; i < l; ++i) {
        index = int(s[i]) - beg;
        if (index < 0 || index >= ALPHA_SIZE) return 0;
    }
    return 1;

}

template <class T, int ALPHA_SIZE>
void trie_type<T,ALPHA_SIZE>::copy(trie_node* source, trie_node* &dest)
{
    if (source == nullptr) return;
    dest = new trie_node();
    dest->end_of_word = source->end_of_word;
    for (int i = 0; i < ALPHA_SIZE; ++i) {
        if (source->children[i] != nullptr) copy(source->children[i],dest->children[i]);
    }
}



template <class T, int ALPHA_SIZE>
void trie_type<T,ALPHA_SIZE>::del(trie_node *node)
{
    if (node == nullptr) return;
    for (int i = 0; i < ALPHA_SIZE; ++i) {
        if (node->children[i] != nullptr) del(node->children[i]);
    }
    delete node;
}


template <class T, int ALPHA_SIZE>
const trie_type<T,ALPHA_SIZE>& trie_type<T,ALPHA_SIZE>::operator =(const trie_type<T,ALPHA_SIZE>& W)
{
    if (head == W.head) return *this;
    if (head != nullptr) del(head);
    copy(W.head,head);
    nwords = W.nwords; beg = W.beg;
    return *this;
}

template <class T, int ALPHA_SIZE>
const trie_type<T,ALPHA_SIZE>& trie_type<T,ALPHA_SIZE>::operator =(trie_type<T,ALPHA_SIZE>&& W)
{
    if (head == W.head) return *this;
    if (head != nullptr) del(head);
    head = W.head; nwords = W.nwords; beg = W.beg;
    W.head = nullptr; W.nwords = 0;
    return *this;
}


template <class T, int ALPHA_SIZE>
void trie_type<T,ALPHA_SIZE>::insert(const T &s)
{
    int i = 0, l = s.length(), index;
    if (!l) return;
    if (!valid_word(s,beg)) return;
    trie_node *current = head;
    if (head == nullptr) {
        head = new trie_node();
        index = int(s[i]) - beg;
        head->children[index] = new trie_node();
        current = head->children[index];
        ++i;
    }

    for (; i < l; ++i) {
        index = int(s[i]) - beg;
        if (current->children[index] == nullptr) {
            current->children[index] = new trie_node();
        }
        current = current->children[index];
    }
    if (current->end_of_word == false) ++nwords;
    current->end_of_word = true;
}




template <class T, int ALPHA_SIZE>
struct trie_type<T,ALPHA_SIZE>::trie_node* trie_type<T,ALPHA_SIZE>::Search(const T& prefix) const
{
    if (!valid_word(prefix,beg) || head == nullptr) return nullptr;
    int l = prefix.length(), index;
    trie_node *current = head;

    for (int i = 0; i < l; ++i) {
        index = int(prefix[i]) - beg;
        if (current->children[index] == nullptr) return nullptr;
        current = current->children[index];
    }

    return current;
}



template <class T, int ALPHA_SIZE>
template <template <typename> class C>
void trie_type<T,ALPHA_SIZE>::appendmatchingword(trie_node *current,C<T> &G,const T &prefix, const T &suffix) const
{
    if (current == nullptr) return;
    if (current->end_of_word) G.push_back(prefix + suffix);

    for (int i = 0; i < ALPHA_SIZE; ++i) {
        if (current->children[i] != nullptr) {
            T temp = suffix;
            temp += char(i + beg);

            appendmatchingword(current->children[i],G,prefix,temp);
        }
    }
}


template <class T, int ALPHA_SIZE>
template <template <typename> class C>
C<T> trie_type<T,ALPHA_SIZE>::getmatchingprefix(const T &prefix) const
{
    C<T> G;
    trie_node *current = Search(prefix);
    appendmatchingword(current,G,prefix);
    return G;
}


template <class T, int ALPHA_SIZE>
bool trie_type<T,ALPHA_SIZE>::containsword(const T& word) const
{
    if (head == nullptr) return false;
    trie_node *node(Search(word));
    if (node == nullptr) return false;
    else return node->end_of_word;
}


} // end namespace val

#endif // TRIE_H_INCLUDED
