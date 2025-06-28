#ifdef __List_h
#define __List_h

/*
====================================================================

LIST.H

0.9	(PeterM, 08-93) Original List and Stack classes in G++.
0.961	(PeterM, JohnV, 08-24-93) Replace Borland list containers with
		vanilla C++ template lists in list.h
0.964	(JohnV, 09-28-93) Add ListIter::delete() and supporting data.

0.97	(JamesA 12-30-93) Add restart() to ListIter
0.971	(JamesA 02-05-94) Add putAt(), getAt(), find() to List, itemName to ListItem.
0.975	(JamesA 03-12-94) Add removeItem() to List.
		(HadonN 06-16-94) Add removeAll(), fix(?) removeItem()
		(HadonN 09-19-94) fix(!) removeItem()
0.976	(Hadon 10-20-94) removed putAt(), getAt(), find() from List
0.977	(Hadon 10-28-94) fix operator=(), (call removeAll())
0.978	(JohnV 06-30-95) add deleteAll() to delete all data and remove all items.
0.979 (JohnV 01-11-96) incorporate undocumented additions (from Acquah?)
1.001 (JohnV 02-29-96) renamed N_items() -> nItems().


1.002 (AndrewL 06-26-25) Mimicking java Arraylist type for efficiency
====================================================================
*/

#include <string.h>

#include <stdlib.h>

template <class Type> class ListItem
{
public:
	ListItem(Type* newitem, ListItem<Type>* nxt, const char* key = "")  // jbha
	{
		data = newitem; next = nxt; itemName = strdup(key);
	}
	~ListItem() { free(itemName); itemName = NULL; next = NULL; data = NULL; }
	/*	~ListItem() { free(itemName); }
		deleteData() {delete data; data = NULL; free(itemName); itemName = NULL;}
		(commented out 01-11-96 JohnV)
	*/
	char* itemName; 					// jbha
	ListItem<Type>* next;
	Type* data;
};


template <class Type> class List
{
public:
	int nitems;
	ListItem<Type>* head;
	List() { head = 0; nitems = 0; }
	List(List<Type>& other) {
		head = other.head; nitems = other.nitems;
		*this = other;
		//MEMORY LEAK
	}
	~List() { removeAll(); }

	void removeItem(Type* anItem);
	void removeAll();
	void deleteAll();
	void add(Type* newitem);

	List<Type>& operator=(List<Type>& other)
	{
		//Not a terribly efficient routine.
		//first copy of list is reversed.
		//reversed_list goes out of scope and is deleted.
		ListIter<Type> walk(other);

		removeAll();
		while (!walk.Done())
		{
			add(walk.current());
			walk++;
		}
		return *this;
	};

	Type* pop(void) {
		if (head)
		{
			Type* retval = head->data;
			ListItem<Type>* to_delete = head;
			head = head->next;
			delete to_delete;
			nitems--;
			return retval;
		}
		else return 0;
	}

	int isEmpty() { return (head == 0); }
	void push(Type* newitem) { add(newitem); };
	int nItems() { return nitems; };
};

template <class Type> void List<Type>::add(Type* newitem)
{
	ListItem<Type>* New = new ListItem<Type>(newitem, head);

	nitems++;
	head = New;
}

template <class Type>
void List<Type>::removeItem(Type* anItem)
// use when you need to find the item first, else use version in ListIter
{
	ListItem<Type>* prevPtr, * currentPtr;

	if (head == NULL) {
		return;
	}
	else if (head->data == anItem) {
		currentPtr = head->next;
		delete head;
		head = currentPtr;
		nitems--;
	}
	else {
		prevPtr = head;
		currentPtr = prevPtr->next;
		while (!(currentPtr == NULL)) {
			if (currentPtr->data == anItem) {
				prevPtr->next = currentPtr->next;
				delete currentPtr;
				currentPtr = prevPtr->next;
				nitems--;
			}
			else {
				prevPtr = currentPtr;
				currentPtr = currentPtr->next;
			}
		}
	}
}

// removeAll: removes all list elements, but NOT the data!

template <class Type>
void List<Type>::removeAll()
{
	ListItem<Type>* next;
	while (head)
	{
		next = head->next;
		delete head;
		head = next;
	}
	nitems = 0;
}

// deleteAll: removes all list elements AND data

template <class Type>
void List<Type>::deleteAll()
{
	ListItem<Type>* next;
	while (head)
	{
		next = head->next;
		delete head->data;
		delete head;
		head = next;
	}
	nitems = 0;
}

/*
	Acquah's version:

  template <class Type>
  void List<Type>::deleteAllData()
  {
	ListItem<Type> *next, *current;
	current = head; 					// necc because deleteData() sometimes
	head = NULL;          				// ... (eventually) calls removeItem()
	while ( current )
	{
		next = (*current).next;
		current->deleteData();
		delete current;
		current = next;
	}
	nitems = 0;
  }

*/

////////////////////////////////////////////////////////////////////////////////

template <class Type> class ListIter
{
public:
	ListItem<Type>* currentPtr;
	ListItem<Type>* previousPtr;
	List<Type>* theList;
	ListIter() {}						// empty constructor, must be used with restart
	ListIter(List<Type>& aList) { theList = &aList; restart(); };
	void	restart() { previousPtr = currentPtr = theList->head; };
	void	restart(List<Type>& aList) { theList = &aList; restart(); };
	void	restart(Type& aGroup);
	Type* data() { if (currentPtr) return currentPtr->data; else return NULL; };
	Type* current() { return data(); };
	Type* operator () () { return data(); };
	void	operator++(int) {
		if (currentPtr) {
			previousPtr = currentPtr;
			currentPtr = currentPtr->next;
		}
		else restart();
	};
	int	Done() { return currentPtr == NULL; };
	int	isEmpty() { return !(theList->head == 0); };
	void	deleteCurrent();
};


//	This function removes an item from the list.  Might be cleaner
//	if it is inside List<Type>, but do not want order n search
//	for element.
template <class Type>
void ListIter<Type>::deleteCurrent()
{
	if (currentPtr == NULL) return;		//	bail out if empty or at end
	//	if at head of list, point head to next item
	if (currentPtr == theList->head)
	{
		theList->head = currentPtr->next;
		previousPtr = NULL;
	}
	else previousPtr->next = currentPtr->next;
	delete	current();						//	delete the data
	delete	currentPtr;						//	delete the list item
	currentPtr = previousPtr;
	theList->nitems--;
}

//  Some useful macros
#define ApplyToList(function,theList,listType) { ListIter<listType> titer(theList); for(titer.restart();!titer.Done();titer++) titer.current()->function;}


#else
#include <string.h>
#include <vector>
#include <stdlib.h>
typedef struct {
	Node* prev;
	Node* next;
	void* data;
} Node;

class List {
public:
	Node* head;
	Node* tail;
	List() 
	{
		head = (Node*)malloc(sizeof(Node));
		tail = head;
	}

	void appendToHead(void* d, size_t data_size) {
		Node* newNode = (Node*)malloc(sizeof(Node));
		newNode->data = malloc(data_size);
		newNode->next = head;

		memcpy(newNode->data, d, data_size);

		head = newNode;
	}

	void appendToTail(void* d, size_t data_size) {
		Node* newNode = (Node*)malloc(sizeof(Node));
		newNode->data = malloc(data_size);
		tail->next = newNode;

		memcpy(newNode->data, d, data_size);

		tail = newNode;
	}
	void nodalRemove() {

	}
};

template <typename T>
class ArrayList : public List 
{
public:
	unsigned int nItems;
	ArrayList() : List() {
		nItems = 0;
	}
	T operator[](unsigned int x) {
		auto curr = head;
		curr += x * sizeof(Node);
		T* result = (T*)(curr->data);
		if (result != nullptr)
			return *result;
		return NULL;
	}
	void add(T* data) 
	{
		appendToTail(data, sizeof(T));
	}
	int remove(T item)
	{
		if (!head || nItems == 0) {
			return -1;
		}
		auto curr = head;
		while (curr && curr->data && *(curr->data) != item) {
			curr = curr->next;
		}
		if (!curr)
		{
			return 1;
		}

		Node* prev = curr->prev;
		prev->next = curr->next;
		free(curr->data);
		free(curr);
		return 0;
	}
	std::vector<T*> removeAll()
	{
		std::vector<T*> survive = {};
		if (!head)
			return survive;
		Node* curr = head;
		while (curr)
		{
			Node* prev = curr;
			curr = curr->next;
			survive.push_back((T*)prev->data);
			free(prev);
		}
		return survive;
	}
	int deleteAll()
	{
		if (!head)
			return -1;
		Node* curr = head;
		while (curr)
		{
			Node* next = curr->next;
			//order matters
			free(curr->data);
			free(curr);
			curr = next;
		}
		head = tail = nullptr;
		nItems = 0;
		return 0;
	}
};

class GenericList : public List 
{
	unsigned int nItems;
	GenericList() : List() {}

	template <typename T>
	void add(T newItem) {
		add((void*)&newItem, sizeof(T));
	}
	void add(void* newItem, size_t size)
	{
		appendToTail(newItem, size);
	}

	template <typename T>
	int remove(T item) {
		auto curr = head;
		while (curr && curr->data && *(curr->data) != item) {
			curr = curr->next;
		}
		if (!curr)
		{
			return 1;
		}

		auto prev = curr->prev;
		prev->next = curr->next;
		free(curr->data);
		free(curr);
		return 0;
	}
	GenericList removeAll()
	{
		GenericList survive;
		if (!head)
			return survive;
		Node* curr = head;
		while (curr)
		{
			Node* prev = curr;
			curr = curr->next;
			survive.add(prev->data);
			free(prev);
		}
		return survive;
	}

	int deleteAll()
	{
		if (!head)
			return -1;
		Node* curr = head;
		while (curr)
		{
			Node* next = curr->next;
			//order matters
			free(curr->data);
			free(curr);
			curr = next;
		}
		head = tail = nullptr;
		nItems = 0;
		return 0;
	}
};

#endif