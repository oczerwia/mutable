#pragma once

#include <cstdint>
#include <iostream>
#include <mutable/mutable-config.hpp>
#include <mutable/parse/AST.hpp>
#include <vector>

// This is for the extract join columns function
#include <map>
#include <set>
#include <string>

namespace m
{

    namespace cnf
    {

        /** A `Predicate` contains a `Expr` of `Boolean` type in either *positive* or *negative* form. */
        struct M_EXPORT Predicate
        {
        private:
            uintptr_t literal_; ///< pointer to Expr; LSB is 1 iff literal is negated

            explicit Predicate(uintptr_t l) : literal_(l) {}

        public:
            /** Creates a *positive* `Predicate` from `e`. */
            static Predicate Positive(const ast::Expr *e) { return Predicate(reinterpret_cast<uintptr_t>(e) | 0x0UL); }
            /** Creates a *negative* `Predicate` from `e`. */
            static Predicate Negative(const ast::Expr *e) { return Predicate(reinterpret_cast<uintptr_t>(e) | 0x1UL); }
            /** Creates a `Predicate` from `e`.  The `Predicate` is *negative* iff `is_negative`. */
            static Predicate Create(const ast::Expr *e, bool is_negative) { return is_negative ? Negative(e) : Positive(e); }

            /** Returns `true` iff this `Predicate` is *negative*. */
            bool negative() const { return literal_ & 0x1UL; }

            /** Returns the `Expr` within this `Predicate`. */
            ast::Expr &expr() { return *reinterpret_cast<ast::Expr *>(literal_ & ~0b11UL); }
            /** Returns the `Expr` within this `Predicate`. */
            const ast::Expr &expr() const { return *reinterpret_cast<const ast::Expr *>(literal_ & ~0b11UL); }
            /** Returns the `Expr` within this `Predicate`. */
            const ast::Expr &operator*() const { return expr(); }
            /** Returns the `Expr` within this `Predicate`. */
            const ast::Expr *operator->() const { return &expr(); }

            /** Returns `true` iff this `Predicate` is an equi-predicate, i.e. an equality comparison of two designators. */
            bool is_equi() const
            {
                auto binary = cast<const ast::BinaryExpr>(&expr());
                if (not binary)
                    return false;
                if (not negative() and binary->tok != TK_EQUAL)
                    return false; // `=` is ok
                if (negative() and binary->tok != TK_BANG_EQUAL)
                    return false; // negated `!=` is ok
                return is<const ast::Designator>(binary->lhs) and is<const ast::Designator>(binary->rhs);
            }
            /** Returns `true` iff this `Predicate` is nullable, i.e. may evaluate to `NULL` at runtime. */
            bool can_be_null() const { return expr().can_be_null(); }

            /** Returns a negated version of this `Predicate`, i.e.\ if this `Predicate` is *positive*, the returned `Predicate`
             * is *negative*. */
            Predicate operator!() const { return Predicate(literal_ ^ 0x1UL); }

            /** Returns `true` iff `other` is equal to `this`.  Two `Predicate`s are equal, iff they have equal `Expr` and the
             * same *sign*. */
            bool operator==(Predicate other) const
            {
                return this->negative() == other.negative() and this->expr() == other.expr();
            }
            /** Returns `true` iff `other` is not equal to `this`.  Two `Predicate`s are equal, iff they have the same
             * `Expr` and the same *sign*. */
            bool operator!=(Predicate other) const { return not operator==(other); }

            /** Compare `Predicate`s by the location of their referenced `Expr` in memory and their sign.  Negative
             * `Predicate`s are larger than positive `Predicate`s of the same expression. */
            bool operator<(Predicate other) const { return this->literal_ < other.literal_; }

            /** Print as SQL expression. */
            void to_sql(std::ostream &out) const;

            /** Print a textual representation of `pred` to `out`. */
            friend std::ostream &M_EXPORT operator<<(std::ostream &out, const Predicate &pred);
            friend M_EXPORT std::string to_string(const Predicate &pred)
            {
                std::ostringstream oss;
                oss << pred;
                return oss.str();
            }

            void dump(std::ostream &out) const;
            void dump() const;

            std::string get_column_name() const
            {
                const auto *designator = cast<const ast::Designator>(&expr());
                if (designator)
                {
                    if (designator->has_table_name())
                    {
                        return std::string(*designator->table_name.text) + "." + std::string(*designator->attr_name.text);
                    }
                    else
                    {
                        return std::string(*designator->attr_name.text);
                    }
                }
                return ""; // Unsupported expression type
            }
        };

        /** A `cnf::Clause` represents a **disjunction** of `Predicate`s. */
        struct M_EXPORT Clause : public std::vector<Predicate>
        {
            using std::vector<Predicate>::vector; // c'tor

            bool operator<=(const Clause &other) const;
            bool operator>=(const Clause &other) const { return other <= *this; }
            bool operator==(const Clause &other) const { return *this >= other and * this <= other; }
            bool operator!=(const Clause &other) const { return not operator==(other); }

            /** Returns a `Schema` instance containing all required definitions (of `Attribute`s and other `Designator`s). */
            Schema get_required() const
            {
                Schema required;
                for (auto &P : *this)
                    required |= P->get_required();
                return required;
            }

            /** Returns `true` iff this `cnf::Clause` formula is an equi-predicate, i.e. a single equality comparison of two
             * designators. */
            bool is_equi() const
            {
                if (size() != 1)
                    return false;
                auto &literal = operator[](0);
                return literal.is_equi();
            }
            /** Returns `true` iff this `cnf::Clause` formula is nullable, i.e. may evaluate to `NULL` at runtime. */
            bool can_be_null() const
            {
                for (auto &P : *this)
                    if (P.can_be_null())
                        return true;
                return false;
            }

            /** Print as SQL expression. */
            void to_sql(std::ostream &out) const;

            /** Print a textual representation of `clause` to `out`. */
            friend std::ostream &M_EXPORT operator<<(std::ostream &out, const Clause &clause);
            friend M_EXPORT std::string to_string(const Clause &clause)
            {
                std::ostringstream oss;
                oss << clause;
                return oss.str();
            }

            void dump(std::ostream &out) const;
            void dump() const;
        };

        /** A `CNF` represents a **conjunction** of `cnf::Clause`s. */
        struct M_EXPORT CNF : public std::vector<Clause>
        {
            using std::vector<Clause>::vector; // c'tor

            /** Returns a `Schema` instance containing all required definitions (of `Attribute`s and other `Designator`s). */
            Schema get_required() const
            {
                Schema required;
                for (auto &clause : *this)
                    required |= clause.get_required();
                return required;
            }

            /** Returns `true` iff this `CNF` formula is an equi-predicate, i.e. conjunction of equality comparisons of each
             * two designators. */
            bool is_equi() const
            {
                if (size() == 0)
                    return false;
                for (auto &clause : *this)
                    if (not clause.is_equi())
                        return false;
                return true;
            }
            /** Returns `true` iff this `CNF` formula is nullable, i.e. may evaluate to `NULL` at runtime. */
            bool can_be_null() const
            {
                for (auto &clause : *this)
                    if (clause.can_be_null())
                        return true;
                return false;
            }

            bool operator<=(const CNF &other) const;
            bool operator>=(const CNF &other) const { return other <= *this; }
            bool operator==(const CNF &other) const { return *this >= other and * this <= other; }
            bool operator!=(const CNF &other) const { return not operator==(other); }

            /** Print as SQL expression. */
            void to_sql(std::ostream &out) const;

            /** Print a textual representation of `cnf` to `out`. */
            friend M_EXPORT std::ostream &operator<<(std::ostream &out, const CNF &cnf);
            friend M_EXPORT std::string to_string(const CNF &cnf)
            {
                std::ostringstream oss;
                oss << cnf;
                return oss.str();
            }

            void dump(std::ostream &out) const;
            void dump() const;

            std::map<std::string, std::set<std::string>> get_join_columns() const
            {
                std::map<std::string, std::set<std::string>> result;
                for (const auto &clause : *this)
                {
                    for (const auto &pred : clause)
                    {
                        if (auto bin = cast<const ast::BinaryExpr>(&pred.expr()))
                        {
                            if (bin->op().type == m::TK_EQUAL)
                            {
                                auto lhs = cast<const ast::Designator>(bin->lhs.get());
                                auto rhs = cast<const ast::Designator>(bin->rhs.get());
                                if (lhs && rhs)
                                {
                                    auto lhs_table_name = std::string(*lhs->table_name.text);
                                    auto lhs_column_name = std::string(*lhs->attr_name.text);

                                    auto rhs_table_name = std::string(*rhs->table_name.text);
                                    auto rhs_column_name = std::string(*rhs->attr_name.text);

                                    result[lhs_table_name].insert(lhs_column_name);
                                    result[rhs_table_name].insert(rhs_column_name);
                                }
                            }
                        }
                    }
                }
                return result;
            }

            std::map<std::string, std::set<std::string>> get_filter_columns() const
            {
                std::map<std::string, std::set<std::string>> result;
                for (const auto &clause : *this)
                {
                    for (const auto &pred : clause)
                    {
                        if (auto bin = cast<const ast::BinaryExpr>(&pred.expr()))
                        {
                            // Extract left side as designator (e.g., table.column)
                            auto lhs = cast<const ast::Designator>(bin->lhs.get());
                            if (lhs && lhs->has_table_name())
                            {
                                auto table_name = std::string(*lhs->table_name.text);
                                auto column_name = std::string(*lhs->attr_name.text);
                                result[table_name].insert(column_name);
                            }

                            // Also check right side in case filter condition has column on right
                            auto rhs = cast<const ast::Designator>(bin->rhs.get());
                            if (rhs && rhs->has_table_name())
                            {
                                auto table_name = std::string(*rhs->table_name.text);
                                auto column_name = std::string(*rhs->attr_name.text);
                                result[table_name].insert(column_name);
                            }
                        }
                    }
                }
                return result;
            }
        };

        /** Returns the **logical or** of two `cnf::Clause`s, i.e.\ the disjunction of the `Predicate`s of `lhs` and `rhs`.
         */
        Clause M_EXPORT operator||(const Clause &lhs, const Clause &rhs);

        /** Returns the **logical and** of two `cnf::Clause`s, i.e.\ a `CNF` with the two `cnf::Clause`s `lhs` and `rhs`.
         */
        CNF M_EXPORT operator&&(const Clause &lhs, const Clause &rhs);

        /** Returns the **logical and** of two `CNF`s, i.e.\ the conjunction of the `cnf::Clause`s of `lhs` and `rhs`. */
        CNF M_EXPORT operator&&(const CNF &lhs, const CNF &rhs);

        /** Returns the **logical or** of two `CNF`s.  It is computed using the [*distributive law* from Boolean
         * algebra](https://en.wikipedia.org/wiki/Distributive_property): *P ∨ (Q ∧ R) ↔ ((P ∨ Q) ∧ (P ∨ R))* */
        CNF M_EXPORT operator||(const CNF &lhs, const CNF &rhs);

        /** Returns the **logical negation** of a `cnf::Clause`.  It is computed using [De Morgan's laws]
         * (https://en.wikipedia.org/wiki/De_Morgan%27s_laws): *¬(P ∨ Q) ↔ ¬P ∧ ¬Q*. */
        CNF M_EXPORT operator!(const Clause &clause);

        /** Returns the **logical negation** of a `CNF`.  It is computed using [De Morgan's laws]
         * (https://en.wikipedia.org/wiki/De_Morgan%27s_laws): *¬(P ∧ Q) ↔ ¬P ∨ ¬Q*. */
        CNF M_EXPORT operator!(const CNF &cnf);

        /** Converts the `Boolean` `Expr` `e` to a `CNF`. */
        CNF M_EXPORT to_CNF(const ast::Expr &e);
        /** Converts the `Boolean` `Expr` of `c` to a `CNF`. */
        CNF M_EXPORT get_CNF(const ast::Clause &c);

    }

}

namespace std
{

    template <>
    struct hash<m::cnf::Predicate>
    {
        std::size_t operator()(m::cnf::Predicate P) const
        {
            return P.expr().hash() ^ (-P.negative());
        }
    };

}
