#ifndef GRAPHENE_LIST
#define GRAPHENE_LIST

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
#include <vector>
#include <stdlib.h>

struct Node {
	Node* prev;
	Node* next;
	void* data;
};

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